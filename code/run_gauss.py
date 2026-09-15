# -*- coding: utf-8 -*-
"""
高斯光束照明下的验证与优化（含优化后的稳定性验证）。

流程
----
1. 按公式重画铺满 1920x1080 的螺旋锥透镜相位；
2. 平顶照明与高斯照明（1/e² 强度半径 w = 4.32 mm）下的 100 mm 光斑基线；
3. 以高斯照明为目标重新优化径向相位（目标：100 mm 处主瓣能量占比最大）；
4. 优化结果验证：高斯照明主结果 + 平顶照明 + 不同束腰 w = 3.0 / 4.32 / 5.5 mm 的稳定性。

输出：gauss_data.npz / gauss_meta.json / phase_optimized_gauss_1920x1080.png
"""

import json
import time
import numpy as np
import slm_common as C

PAD = 64
Z = 0.100
R_C = C.R_ZERO
W_MAIN = 2.5e-3                  # 1/e² 强度半径，对应光束直径 5 mm
W_LIST = [1.8e-3, 2.5e-3, 3.5e-3]
NODE_STEP = 8e-6
ITERS = 200
LR = 0.05

t0 = time.time()
phi0, res_rms, phi_in = C.fit_phi0()
X, Y = C.grid_coords()
R = np.hypot(X, Y)
TH = np.arctan2(Y, X)
BASE = -C.K * C.SIN_A * R
phase_full = np.mod(BASE + C.ELL * TH + phi0, C.FULL)

nodes = np.arange(0.0, C.R_CORNER + NODE_STEP, NODE_STEP)
NN = nodes.size
tt = np.clip(R / NODE_STEP, 0, NN - 1.001)
ni = np.floor(tt).astype(int)
fr = tt - ni

HH, WW = C.SLM_H + 2 * PAD, C.SLM_W + 2 * PAD
H = C.transfer_function(Z, WW, HH)
HI = np.conj(H)
yy, xx = np.mgrid[0:HH, 0:WW].astype(float)
CEN = (C.SLM_W / 2.0 + PAD - 0.5, C.SLM_H / 2.0 + PAD - 0.5)
RR = np.hypot((xx - CEN[0]) * C.DX, (yy - CEN[1]) * C.DX)
CHI = (RR <= R_C).astype(float)


def fwd(E):
    Ep = np.zeros((HH, WW), np.complex128)
    Ep[PAD:PAD + C.SLM_H, PAD:PAD + C.SLM_W] = E
    return np.fft.ifft2(np.fft.fft2(Ep) * H)


def bwd(V):
    U = np.fft.ifft2(np.fft.fft2(V) * HI)
    return U[PAD:PAD + C.SLM_H, PAD:PAD + C.SLM_W]


def field(delta, gauss_w):
    phi = BASE + C.ELL * TH + delta[ni] * (1 - fr) + delta[ni + 1] * fr
    return C.incident_amplitude(X, Y, gauss_w) * np.exp(1j * phi)


def metrics(I):
    cm = C.core_metrics(I, CEN)
    _, main_r, rings = C.sidelobe_table(I, CEN)
    ee = C.encircled(I, CEN, (20e-6, 37.2e-6, 50e-6, 100e-6, 200e-6))
    return dict(peak=float(I.max()), r_peak=float(main_r * 1e6),
                r_zero=float(cm["r_zero"] * 1e6), eta_main=float(cm["eta_main"]),
                rings=[float(db) for _, _, _, _, db in rings],
                ee={f"{int(q*1e6)}": float(v) for q, v in ee.items()})


print(f"[1] 基线（平顶与高斯照明）  t={time.time()-t0:.0f}s")
base = {}
I_flat = np.abs(fwd(field(np.zeros(NN), None))) ** 2
base["flat"] = metrics(I_flat)
print(f"    平顶照明        峰值 {base['flat']['peak']:9.1f}  环半径 "
      f"{base['flat']['r_peak']:5.1f} µm  第一旁瓣 {base['flat']['rings'][0]:6.2f} dB  "
      f"主瓣能量 {base['flat']['eta_main']*100:6.3f}%")
for tag, w in (("w1.8", 1.8e-3), ("w2.5", 2.5e-3), ("w3.5", 3.5e-3)):
    I = np.abs(fwd(field(np.zeros(NN), w))) ** 2
    base[tag] = metrics(I)
    print(f"    高斯 w={w*1e3:.2f} mm  峰值 {base[tag]['peak']:9.1f}  环半径 "
          f"{base[tag]['r_peak']:5.1f} µm  第一旁瓣 {base[tag]['rings'][0]:6.2f} dB  "
          f"主瓣能量 {base[tag]['eta_main']*100:6.3f}%", flush=True)


def obj_grad(delta, want_grad=True):
    E = field(delta, W_MAIN)
    U = fwd(E)
    I = np.abs(U) ** 2
    D = float(I.sum())
    J = float((CHI * I).sum() / D)
    if not want_grad:
        return J, None
    g = bwd(CHI * U / D)
    gp = 2.0 * np.imag(np.conj(E) * g)
    gn = np.bincount(ni.ravel(), weights=(gp * (1 - fr)).ravel(), minlength=NN)
    gn += np.bincount((ni + 1).ravel(), weights=(gp * fr).ravel(), minlength=NN)
    return J, gn


print(f"[2] 高斯照明 (w = {W_MAIN*1e3:.2f} mm) 下的相位优化  t={time.time()-t0:.0f}s")
delta = np.zeros(NN)
J0, g0 = obj_grad(delta)
rng = np.random.default_rng(0)
v = rng.normal(size=NN)
v /= np.linalg.norm(v)
hs = 2e-4
fd = (obj_grad(hs * v, False)[0] - obj_grad(-hs * v, False)[0]) / (2 * hs)
an = float(np.dot(g0, v))
print(f"    初始主瓣能量 {J0*100:.3f}%；伴随梯度与有限差分偏差 "
      f"{abs(fd-an)/abs(fd)*100:.3f}%")

m_ad, v_ad = np.zeros(NN), np.zeros(NN)
b1, b2, eps = 0.9, 0.999, 1e-8
hist, best_J, best_d = [], J0, delta.copy()
for it in range(1, ITERS + 1):
    J, g = obj_grad(delta)
    m_ad = b1 * m_ad + (1 - b1) * g
    v_ad = b2 * v_ad + (1 - b2) * g ** 2
    delta = delta + LR * (m_ad / (1 - b1 ** it)) / (np.sqrt(v_ad / (1 - b2 ** it)) + eps)
    hist.append(float(J))
    if J > best_J:
        best_J, best_d = float(J), delta.copy()
    if it % 20 == 0 or it == 1:
        print(f"    it {it:4d}  主瓣能量 {J*100:6.3f}%  t={time.time()-t0:.0f}s", flush=True)
print(f"    优化结束：主瓣能量 {hist[0]*100:.3f}% → {best_J*100:.3f}%")

print(f"[3] 优化结果验证  t={time.time()-t0:.0f}s")
I_opt = np.abs(fwd(field(best_d, W_MAIN))) ** 2
opt_main = metrics(I_opt)
verify = {}
for tag, w in (("flat", None), ("w1.8", 1.8e-3), ("w2.5", 2.5e-3), ("w3.5", 3.5e-3)):
    I = np.abs(fwd(field(best_d, w))) ** 2
    verify[tag] = metrics(I)
    print(f"    优化相位 + {tag:5s}  峰值 {verify[tag]['peak']:9.1f}  环半径 "
          f"{verify[tag]['r_peak']:5.1f} µm  第一旁瓣 {verify[tag]['rings'][0]:6.2f} dB  "
          f"主瓣能量 {verify[tag]['eta_main']*100:6.2f}%", flush=True)

phi_opt = np.mod(BASE + C.ELL * TH + best_d[ni] * (1 - fr) + best_d[ni + 1] * fr, C.FULL)
C.phase_to_png(phi_opt, "phase_optimized_gauss_1920x1080.png")
C.phase_to_png(phase_full, "phase_full_1920x1080.png")

np.savez_compressed(
    "gauss_data.npz",
    phase_full=phase_full.astype(np.float32), phi_opt=phi_opt.astype(np.float32),
    phi_in=phi_in.astype(np.float32), delta=best_d, hist=np.array(hist),
    R=R.astype(np.float32), R_C=R_C, CEN=np.array(CEN),
    I_flat=I_flat.astype(np.float32), I_opt=I_opt.astype(np.float32),
)
json.dump(dict(phi0=float(phi0), phi0_res_rms=float(res_rms),
               alpha_deg=float(np.degrees(C.ALPHA)), sin_a=float(C.SIN_A),
               r_peak_um=float(C.R_PEAK*1e6), r_zero_um=float(C.R_ZERO*1e6),
               r_ring_um=float(C.R_RING*1e6), z_mm=Z*1e3,
               gauss_w_um=W_MAIN*1e6, iters=ITERS, lr=LR,
               J_init=float(hist[0]), J_best=float(best_J),
               base=base, verify=verify, opt=opt_main,
               node_count=int(NN), delta_min=float(best_d.min()),
               delta_max=float(best_d.max())),
          open("gauss_meta.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"[3] 完成，用时 {time.time()-t0:.0f}s")
