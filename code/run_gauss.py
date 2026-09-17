# -*- coding: utf-8 -*-
"""
高斯光束照明下的验证与优化（含优化后的稳定性验证）。

流程
----
1. 按公式重画铺满 1920x1080 的螺旋锥透镜相位（见 optimize.RadialPhaseOptimizer）；
2. 平顶照明与高斯照明（1/e² 强度直径 5 mm，w = 2.5 mm）下的 100 mm 光斑基线；
3. 以高斯照明为目标重新优化径向相位（调用 optimize 模块，目标：100 mm 处主瓣能量占比最大）；
4. 优化结果验证：高斯照明主结果 + 平顶照明 + 不同光束直径的稳定性。

输出：gauss_data.npz / gauss_meta.json / phase_optimized_gauss_1920x1080.png
"""

import json
import time
import numpy as np
import slm_common as C
import optimize as O

W_MAIN = O.W_MAIN
W_LIST = [1.8e-3, 2.5e-3, 3.5e-3]
Z = O.Z
R_C = O.R_C

t0 = time.time()
opt = O.RadialPhaseOptimizer()
phase_full = opt.phase_full


def metrics(I):
    cm = C.core_metrics(I, opt.CEN)
    _, main_r, rings = C.sidelobe_table(I, opt.CEN)
    ee = C.encircled(I, opt.CEN, (20e-6, 37.2e-6, 50e-6, 100e-6, 200e-6))
    return dict(peak=float(I.max()), r_peak=float(main_r * 1e6),
                r_zero=float(cm["r_zero"] * 1e6), eta_main=float(cm["eta_main"]),
                rings=[float(d) for _, _, _, _, d in rings],
                ee={f"{int(q*1e6)}": float(v) for q, v in ee.items()})


# ---------------------------------------------------------------- [1] 基线
print(f"[1] 基线（平顶与高斯照明）  t={time.time()-t0:.0f}s", flush=True)
base = {}
I_flat = np.abs(opt.fwd(opt.field(np.zeros(opt.NN), None))) ** 2
base["flat"] = metrics(I_flat)
print(f"    平顶照明        峰值 {base['flat']['peak']:9.1f}  环半径 "
      f"{base['flat']['r_peak']:5.1f} µm  第一旁瓣 {base['flat']['rings'][0]:6.2f} dB  "
      f"主瓣能量 {base['flat']['eta_main']*100:6.3f}%", flush=True)
for tag, w in (("w1.8", 1.8e-3), ("w2.5", 2.5e-3), ("w3.5", 3.5e-3)):
    I = np.abs(opt.fwd(opt.field(np.zeros(opt.NN), w))) ** 2
    base[tag] = metrics(I)
    print(f"    高斯 w={w*1e3:.2f} mm  峰值 {base[tag]['peak']:9.1f}  环半径 "
          f"{base[tag]['r_peak']:5.1f} µm  第一旁瓣 {base[tag]['rings'][0]:6.2f} dB  "
          f"主瓣能量 {base[tag]['eta_main']*100:6.3f}%", flush=True)

# ---------------------------------------------------------------- [2] 优化
print(f"[2] 高斯照明（w = {W_MAIN*1e3:.2f} mm）下的相位优化  t={time.time()-t0:.0f}s",
      flush=True)
best_d, hist, best_J, chk = opt.run()
print(f"    优化结束：主瓣能量 {chk['J0']*100:.3f}% → {best_J*100:.3f}%", flush=True)

# ---------------------------------------------------------------- [3] 验证
print(f"[3] 优化结果验证  t={time.time()-t0:.0f}s", flush=True)
I_opt = np.abs(opt.fwd(opt.field(best_d, W_MAIN))) ** 2
opt_main = metrics(I_opt)
verify = {}
for tag, w in (("flat", None), ("w1.8", 1.8e-3), ("w2.5", 2.5e-3), ("w3.5", 3.5e-3)):
    I = np.abs(opt.fwd(opt.field(best_d, w))) ** 2
    verify[tag] = metrics(I)
    print(f"    优化相位 + {tag:5s}  峰值 {verify[tag]['peak']:9.1f}  环半径 "
          f"{verify[tag]['r_peak']:5.1f} µm  第一旁瓣 {verify[tag]['rings'][0]:6.2f} dB  "
          f"主瓣能量 {verify[tag]['eta_main']*100:6.2f}%", flush=True)

# ---------------------------------------------------------------- [4] 保存
phi_opt = np.mod(opt.BASE + C.ELL * opt.TH
                 + best_d[opt.ni] * (1 - opt.fr) + best_d[opt.ni + 1] * opt.fr, C.FULL)
C.phase_to_png(phi_opt, "phase_optimized_gauss_1920x1080.png")
C.phase_to_png(phase_full, "phase_full_1920x1080.png")

np.savez_compressed(
    "gauss_data.npz",
    phase_full=phase_full.astype(np.float32), phi_opt=phi_opt.astype(np.float32),
    phi_in=opt.phi_in.astype(np.float32), delta=best_d, hist=hist,
    R=opt.R.astype(np.float32), R_C=R_C, CEN=np.array(opt.CEN),
    I_flat=I_flat.astype(np.float32), I_opt=I_opt.astype(np.float32),
)
json.dump(dict(phi0=float(opt.phi0), phi0_res_rms=float(opt.res_rms),
               alpha_deg=float(np.degrees(C.ALPHA)), sin_a=float(C.SIN_A),
               r_peak_um=float(C.R_PEAK * 1e6), r_zero_um=float(C.R_ZERO * 1e6),
               r_ring_um=float(C.R_RING * 1e6), z_mm=Z * 1e3,
               gauss_w_um=W_MAIN * 1e6, iters=O.ITERS, lr=O.LR,
               J_init=float(chk["J0"]), J_best=float(best_J),
               grad_err_pct=float(chk["rel_err_pct"]),
               base=base, verify=verify, opt=opt_main,
               node_count=int(opt.NN), delta_min=float(best_d.min()),
               delta_max=float(best_d.max())),
          open("gauss_meta.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"[4] 完成，用时 {time.time()-t0:.0f}s")
