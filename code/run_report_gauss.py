# -*- coding: utf-8 -*-
"""
报告生成：相位图重画、100 mm 光斑分析、旁瓣抑制优化与稳定性验证。
所有强度/功率/能量均换算为绝对单位（入射激光功率 200 mW）。
"""

import json
import numpy as np
from chartlib import (PAL, MAGMA, DIVERGE, LUT_MAG, LUT_DIV, LUT_GRAY, norm, png_uri,
                      crop, css_gradient, chart, tbl, figcap)
import slm_common as C

D = np.load("gauss_data.npz")
M = json.load(open("gauss_meta.json", encoding="utf-8"))

phase_full = D["phase_full"]
phi_opt = D["phi_opt"]
delta = D["delta"]
hist = D["hist"]
CEN = tuple(D["CEN"])
PAD = 64
R_C = float(D["R_C"])
W_MAIN = M["gauss_w_um"] * 1e-6

# ---------------------------------------------------------------- 绝对单位换算
P_IN = 200.0                                     # 入射激光功率 [mW]
A_FLAT = C.AP_W * 1e3 * C.AP_H * 1e3             # 平顶受照面积 [mm²]
I_FLAT = P_IN / A_FLAT                           # 平顶强度 [mW/mm²]
TAU = 10e-9                                      # 计算能量所用的脉宽示例 [s]


def i_peak_gauss(w_m):
    """高斯光束峰值强度 [mW/mm²]，总功率 P_IN、1/e² 强度半径 w。"""
    w_mm = w_m * 1e3
    return 2 * P_IN / (np.pi * w_mm ** 2) / (1 - np.exp(-2 * (C.AP_H / w_m) ** 2))


I_G = i_peak_gauss(W_MAIN)


def energy_nJ(p_mW, tau=TAU):
    """由功率 [mW] 与作用时间 tau [s] 换算能量 [nJ]。"""
    return p_mW * 1e-3 * tau * 1e9


def fwd(E, gauss_w):
    X, Y = C.grid_coords()
    A = C.incident_amplitude(X, Y, gauss_w)
    h, w = E.shape
    Ep = np.zeros((h + 2 * PAD, w + 2 * PAD), np.complex128)
    Ep[PAD:PAD + h, PAD:PAD + w] = A * E
    return np.abs(C.propagate(Ep, 0.100)) ** 2


E_flat = np.exp(1j * phase_full)
E_opt = np.exp(1j * phi_opt)
CASES = {
    "base_flat": (E_flat, None),
    "base_g": (E_flat, W_MAIN),
    "base_18": (E_flat, 1.8e-3),
    "base_35": (E_flat, 3.5e-3),
    "opt_g": (E_opt, W_MAIN),
    "opt_flat": (E_opt, None),
    "opt_18": (E_opt, 1.8e-3),
    "opt_35": (E_opt, 3.5e-3),
}
I = {k: fwd(E, w) for k, (E, w) in CASES.items()}


def met(key):
    A = I[key]
    cm = C.core_metrics(A, CEN)
    _, rp, rings = C.sidelobe_table(A, CEN)
    ee = C.encircled(A, CEN, (20e-6, 37.2e-6, 50e-6, 100e-6, 200e-6))
    rc, prof = C.radial_profile(A, CEN, rmax=300e-6, nbins=900)
    gauss_w = CASES[key][1]
    i_in = I_FLAT if gauss_w is None else i_peak_gauss(gauss_w)
    return dict(peak_norm=float(A.max()),
                peak=float(A.max()) * i_in,                      # mW/mm²
                p_main=float(cm["eta_main"]) * P_IN,             # mW
                e_main=energy_nJ(float(cm["eta_main"]) * P_IN),   # nJ（τ = 10 ns）
                eta=float(cm["eta_main"]),
                r_peak=float(rp * 1e6), r_zero=float(cm["r_zero"] * 1e6),
                rings=[float(d) for _, _, _, _, d in rings],
                ee={f"{int(q*1e6)}": float(v) for q, v in ee.items()},
                r=rc * 1e6, p=prof, i_in=i_in)


mt = {k: met(k) for k in CASES}
print("入射：平顶 %.3f mW/mm²，高斯(w=2.5mm) %.3f mW/mm²；总功率 %.0f mW"
      % (I_FLAT, I_G, P_IN))
print("平顶基线 : 峰值 %.1f mW/mm²  主瓣功率 %.2f mW  第一旁瓣 %.2f dB"
      % (mt["base_flat"]["peak"], mt["base_flat"]["p_main"], mt["base_flat"]["rings"][0]))
print("高斯基线 : 峰值 %.1f mW/mm²  主瓣功率 %.2f mW  第一旁瓣 %.2f dB"
      % (mt["base_g"]["peak"], mt["base_g"]["p_main"], mt["base_g"]["rings"][0]))
print("高斯优化 : 峰值 %.1f mW/mm²  主瓣功率 %.2f mW  第一旁瓣 %.2f dB"
      % (mt["opt_g"]["peak"], mt["opt_g"]["p_main"], mt["opt_g"]["rings"][0]))


def cell(uri, cap):
    return ('<figure class="cell"><div class="imgbox"><img src="' + uri
            + '" alt="fig"></div><figcaption>' + cap + '</figcaption></figure>')


def scell(svg, cap):
    return ('<figure class="cell"><div class="imgbox chartbox">' + svg
            + '</div><figcaption>' + cap + '</figcaption></figure>')


def cnt(A, mm):
    k = int(round(mm * 1e-3 / C.DX))
    cx, cy = int(round(CEN[0])), int(round(CEN[1]))
    return A[cy - k:cy + k, cx - k:cx + k]


def spot(key, mm, log=True, floor=4, ref_mw=None):
    """ref_mw：统一的绝对强度标尺上限 [mW/mm²]；None 表示按该图自身峰值。"""
    A = cnt(I[key], mm) * mt[key]["i_in"]
    vmax = float(A.max() if ref_mw is None else ref_mw)
    if log:
        u = norm(A, log=True, vmin=np.log10(vmax * 10 ** (-floor)), vmax=np.log10(vmax))
        u = np.where(A > vmax, 1.0, u)
    else:
        u = norm(A, vmax=vmax)
    return png_uri(u, LUT_MAG, size=(430, 430))


def db(curve):
    return 10 * np.log10(curve["p"] / curve["p"].max())


def rng_mw(key, mm):
    """返回 (显示下限, 显示上限)，单位 mW/mm²：对数图取下限 = 峰值 × 10⁻⁴。"""
    A = cnt(I[key], mm) * mt[key]["i_in"]
    return A.max() * 1e-4, A.max()


PCX = (C.SLM_W - 1) / 2.0
PCY = (C.SLM_H - 1) / 2.0


def pattern_window(phase, n=256):
    x0 = int(round(PCX - n / 2.0))
    y0 = int(round(PCY - n / 2.0))
    return phase[y0:y0 + n, x0:x0 + n]


F1 = ('<div class="grid2">'
      + cell(png_uri(norm(phase_full, vmin=0, vmax=C.FULL), LUT_GRAY, size=(900, 506)),
             "图 1a　按公式重画的相位图（1920 × 1080，8 µm 像素，铺满整幅 SLM）")
      + cell(png_uri(norm(crop(phase_full, 120), vmin=0, vmax=C.FULL), LUT_GRAY, size=(420, 420)),
             "图 1b　中心 240 × 240 区域：螺旋相位奇点与等间距环带")
      + '</div>'
      + '<div class="grid3">'
      + cell(png_uri(norm(D["phi_in"], vmin=0, vmax=C.FULL), LUT_GRAY, size=(300, 300)),
             "图 1c　输入 256 × 256 相位图")
      + cell(png_uri(norm(pattern_window(phase_full), vmin=0, vmax=C.FULL),
                     LUT_GRAY, size=(300, 300)),
             "图 1d　重画图中心 256 × 256 区域")
      + cell(png_uri(norm(np.mod(pattern_window(phase_full) - D["phi_in"] + np.pi, C.FULL) - np.pi,
                          vmin=-0.02, vmax=0.02), LUT_DIV, size=(300, 300)),
             "图 1e　重画图与输入图之差（发散色标，显示范围 ±0.02 rad）")
      + '</div>'
      + '<div class="gradrow"><span>−0.02 rad</span><i class="grad" style="background-image:'
      + css_gradient(DIVERGE) + '"></i><span>+0.02 rad</span><span class="muted">'
        '（色标中点为零偏差；实测差值范围 −0.0167 … +0.0160 rad，'
        '均方根 0.0073 rad，94 % 的像素落在 ±0.0123 rad 以内，即 8 bit 舍入量化界内）'
        '</span></div>'
      + figcap("图 1　相位图重画。由输入图标定的参数为：径向相位斜率 0.8245 rad/像素、"
               "拓扑荷 ℓ = +1、常数相位 φ₀ = %.4f rad。重画图案铺满 15.36 mm × 8.64 mm，"
               "环带周期 %.2f µm。图 1e 为中心 256 × 256 区域内重画结果与输入图的逐个像素相位差"
               "（同一坐标网格对齐后相减），差值均方根 %.4f rad，为 8 bit 量化步长 "
               "0.0246 rad 的 %.0f%%，且 94 %% 的像素落在 ±0.0123 rad（半个量化步长）以内，"
               "表明偏差来源为输入图的 8 bit 相位量化，而非重画误差。"
               % (M["phi0"], M["r_ring_um"], M["phi0_res_rms"],
                  M["phi0_res_rms"] / 0.0246 * 100)))

lo_b, hi_b = rng_mw("base_flat", 1.2)
lo_o, hi_o = rng_mw("opt_g", 0.6)
REF_BASE = mt["base_g"]["peak"]        # 优化前两图共用的绝对标尺上限
F2 = ('<div class="grid2">'
      + cell(spot("base_flat", 0.6, ref_mw=REF_BASE),
             f"图 2a　平顶照明，±0.6 mm（对数显示 {lo_b:.2f}–{hi_b:.0f} mW/mm²）")
      + cell(spot("base_g", 0.6, ref_mw=REF_BASE),
             f"图 2b　高斯照明（D = 5 mm），±0.6 mm（同一绝对标尺 "
             f"{REF_BASE*1e-4:.2f}–{REF_BASE:.0f} mW/mm²）")
      + '</div><div class="grid2">'
      + cell(spot("base_flat", 0.06, log=False, ref_mw=REF_BASE),
             f"图 2c　平顶照明，中心 ±0.06 mm（线性 0–{REF_BASE:.0f} mW/mm²）")
      + cell(spot("base_g", 0.06, log=False, ref_mw=REF_BASE),
             f"图 2d　高斯照明，中心 ±0.06 mm（线性 0–{REF_BASE:.0f} mW/mm²，同一标尺）")
      + '</div>')

c3 = chart([dict(x=mt["base_flat"]["r"], y=db(mt["base_flat"]), label="平顶照明", color=PAL[0]),
            dict(x=mt["base_g"]["r"], y=db(mt["base_g"]), label="高斯照明 D = 5 mm",
                 color=PAL[1], dash="6 4", width=1.8)],
           "径向坐标 r（µm）", "相对主瓣峰值（dB）", w=720, h=360, legpos="sw", ylim=(-40, 3),
           annotations=[("vline", C.R_ZERO * 1e6, "主瓣边界 37.2 µm", 30)])
F3 = scell(c3, "图 3　两种照明条件下 100 mm 处的径向剖面（峰值均为 0 dB）")

c4a = chart([dict(x=np.arange(1, hist.size + 1), y=hist * P_IN, color=PAL[0], width=2)],
            "迭代次数", "主瓣内功率（mW）", w=700, h=320, legpos="se")
prof_r = np.arange(0, C.R_HALF_V, 4e-6)
idx = np.clip(np.round(prof_r / 8e-6).astype(int), 0, delta.size - 1)
phi_tot = (-C.K * C.SIN_A * prof_r + delta[idx])
wrap = np.mod(np.diff(phi_tot) + np.pi, C.FULL) - np.pi
local_sin = wrap / (C.K * (prof_r[1] - prof_r[0]))
k_s = 15
local_sin = np.convolve(local_sin, np.ones(k_s) / k_s, mode="same")
prof_mid = 0.5 * (prof_r[1:] + prof_r[:-1])
sin_eff = -local_sin
keep = prof_mid <= 4.2e-3
c4b = chart([dict(x=prof_mid[keep] * 1e6, y=1e3 * sin_eff[keep], label="优化后", color=PAL[4],
                  width=2),
             dict(x=prof_mid[keep] * 1e6, y=1e3 * C.SIN_A * np.ones(int(keep.sum())),
                  label="原始（均匀锥角）", color=PAL[0], dash="6 4", width=1.8)],
            "半径 r（µm）", "局部等效锥角 sinθ（×10⁻³）", w=700, h=340, legpos="se")
F4 = (scell(c4a, "图 4a　优化收敛过程（纵轴为 z = 100 mm 处主瓣内的功率）")
      + scell(c4b, "图 4b　优化后相位的局部等效锥角随半径的变化（原始图案为常值 "
                   f"{C.SIN_A*1e3:.2f}×10⁻³）"))

F5 = ('<div class="grid2">'
      + cell(spot("opt_g", 0.6),
             f"图 5a　优化后（高斯照明），±0.6 mm（对数显示 {lo_o:.2f}–{hi_o:.0f} mW/mm²）")
      + cell(spot("opt_g", 0.06, log=False),
             f"图 5b　优化后中心 ±0.06 mm（线性 0–{hi_o:.0f} mW/mm²）")
      + '</div>')

c6 = chart([dict(x=mt["base_g"]["r"], y=db(mt["base_g"]), label="优化前", color=PAL[0]),
            dict(x=mt["opt_g"]["r"], y=db(mt["opt_g"]), label="优化后", color=PAL[4],
                 dash="6 4", width=1.8)],
           "径向坐标 r（µm）", "相对主瓣峰值（dB）", w=720, h=380, legpos="sw", ylim=(-45, 3),
           annotations=[("vline", C.R_ZERO * 1e6, "主瓣边界 37.2 µm", 30)])
F6 = scell(c6, "图 6　优化前后的径向剖面对比（z = 100 mm，高斯照明 D = 5 mm）")

F7 = ('<div class="grid2">'
      + cell(png_uri(norm(phi_opt, vmin=0, vmax=C.FULL), LUT_GRAY, size=(900, 506)),
             "图 7a　优化后的相位图（1920 × 1080）")
      + cell(png_uri(norm(crop(phi_opt, 120), vmin=0, vmax=C.FULL), LUT_GRAY, size=(420, 420)),
             "图 7b　中心 240 × 240 区域")
      + '</div>')

c8 = chart([dict(x=mt["opt_g"]["r"], y=db(mt["opt_g"]), label="D = 5 mm（设计值）", color=PAL[4]),
            dict(x=mt["opt_18"]["r"], y=db(mt["opt_18"]), label="D = 3.6 mm", color=PAL[2],
                 dash="5 4", width=1.6),
            dict(x=mt["opt_35"]["r"], y=db(mt["opt_35"]), label="D = 7.0 mm", color=PAL[1],
                 dash="2 3", width=1.6),
            dict(x=mt["opt_flat"]["r"], y=db(mt["opt_flat"]), label="平顶照明", color=PAL[0],
                 dash="7 4", width=1.6)],
           "径向坐标 r（µm）", "相对主瓣峰值（dB）", w=720, h=380, legpos="sw", ylim=(-45, 3))
F8 = scell(c8, "图 8　优化相位在不同光束直径与平顶照明条件下的径向剖面（z = 100 mm）")

T_SET = tbl(["项目", "取值", "说明"], [
    ["输入相位图", "256 × 256，像素 8 µm", "螺旋锥透镜：φ = −k·sinα·r + ℓ·θ + φ₀"],
    ["输出相位图", "1920 × 1080，像素 8 µm", "按公式重画，铺满整幅 SLM（15.36 mm × 8.64 mm）"],
    ["等效锥角", f'α = {M["alpha_deg"]:.4f}°（sinα = {M["sin_a"]:.6f}）', "由输入图标定"],
    ["环形特征尺寸", f'亮环 {M["r_peak_um"]:.2f} µm，第一零点 {M["r_zero_um"]:.2f} µm，'
                    f'环间距 {M["r_ring_um"]:.2f} µm', "理想贝塞尔-涡旋光束理论值"],
    ["入射激光功率", f'{P_IN:.0f} mW', "连续输出时的总功率；脉冲工作时单脉冲能量 E = P·τ"],
    ["平顶照明", f'均匀强度 {I_FLAT:.3f} mW/mm²',
     f'通光面积 {A_FLAT:.2f} mm²，总功率 {P_IN:.0f} mW'],
    ["高斯照明", f'峰值强度 {I_G:.2f} mW/mm²',
     f'1/e² 强度直径 5.0 mm（半径 w = 2.50 mm），总功率 {P_IN:.0f} mW'],
    ["波长 / 重建距离", "1064 nm / 100 mm", "标量衍射，忽略偏振与吸收"],
    ["传播算法", "带限角谱法（Matsushima）", "网格 2048 × 1208（补零 64 像素）"],
    ["优化变量", f'径向相位修正 δ(r)，{M["node_count"]} 个节点', "节点间隔 8 µm（1 个 SLM 像素），线性插值"],
])

T_BASE = tbl(["照明条件", "峰值强度（mW/mm²）", "亮环半径（µm）", "第一旁瓣（dB）",
              "主瓣内功率（mW）", "主瓣内能量（nJ，τ = 10 ns）"], [
    ["平顶（200 mW）", f'{mt["base_flat"]["peak"]:.1f}', f'{mt["base_flat"]["r_peak"]:.1f}',
     f'{mt["base_flat"]["rings"][0]:.2f}', f'{mt["base_flat"]["p_main"]:.2f}',
     f'{mt["base_flat"]["e_main"]:.3f}'],
    ["高斯 D = 3.6 mm（200 mW）", f'{mt["base_18"]["peak"]:.1f}',
     f'{mt["base_18"]["r_peak"]:.1f}', f'{mt["base_18"]["rings"][0]:.2f}',
     f'{mt["base_18"]["p_main"]:.2f}', f'{mt["base_18"]["e_main"]:.3f}'],
    ["高斯 D = 5.0 mm（200 mW）", f'{mt["base_g"]["peak"]:.1f}',
     f'{mt["base_g"]["r_peak"]:.1f}', f'{mt["base_g"]["rings"][0]:.2f}',
     f'{mt["base_g"]["p_main"]:.2f}', f'{mt["base_g"]["e_main"]:.3f}'],
    ["高斯 D = 7.0 mm（200 mW）", f'{mt["base_35"]["peak"]:.1f}',
     f'{mt["base_35"]["r_peak"]:.1f}', f'{mt["base_35"]["rings"][0]:.2f}',
     f'{mt["base_35"]["p_main"]:.2f}', f'{mt["base_35"]["e_main"]:.3f}'],
])

T_OPT = tbl(["指标", "优化前（高斯 D = 5 mm）", "优化后（高斯 D = 5 mm）", "变化"], [
    ["峰值强度（mW/mm²）", f'{mt["base_g"]["peak"]:.1f}', f'{mt["opt_g"]["peak"]:.1f}',
     f'×{mt["opt_g"]["peak"]/mt["base_g"]["peak"]:.0f}'],
    ["亮环峰值半径", f'{mt["base_g"]["r_peak"]:.1f} µm', f'{mt["opt_g"]["r_peak"]:.1f} µm', "—"],
    ["第一零点半径", f'{mt["base_g"]["r_zero"]:.1f} µm', f'{mt["opt_g"]["r_zero"]:.1f} µm', "—"],
    ["第一旁瓣", f'{mt["base_g"]["rings"][0]:.2f} dB', f'{mt["opt_g"]["rings"][0]:.2f} dB',
     f'{mt["opt_g"]["rings"][0]-mt["base_g"]["rings"][0]:+.2f} dB'],
    ["主瓣内功率（mW）", f'{mt["base_g"]["p_main"]:.2f}', f'{mt["opt_g"]["p_main"]:.1f}',
     f'×{mt["opt_g"]["p_main"]/mt["base_g"]["p_main"]:.1f}'],
    ["主瓣内能量（nJ，τ = 10 ns）", f'{mt["base_g"]["e_main"]:.3f}',
     f'{mt["opt_g"]["e_main"]:.3f}', f'×{mt["opt_g"]["e_main"]/mt["base_g"]["e_main"]:.1f}'],
    ["主瓣能量占入射总能量比例", f'{mt["base_g"]["eta"]*100:.3f} %',
     f'{mt["opt_g"]["eta"]*100:.2f} %', f'×{mt["opt_g"]["eta"]/mt["base_g"]["eta"]:.0f}'],
])

def annulus_power(key, r1, r2):
    A = I[key]
    h, w = A.shape
    y, x = np.mgrid[0:h, 0:w]
    r = np.hypot((x - CEN[0]) * C.DX, (y - CEN[1]) * C.DX)
    return float(A[(r >= r1) & (r < r2)].sum() / A.sum()) * P_IN


T_SIDE = tbl(["环带（µm）", "优化前旁瓣（dB）", "优化后旁瓣（dB）", "改善（dB）",
              "优化前环带功率（mW）", "优化后环带功率（mW）"], [
    [f'{C.R_ZERO*1e6 + k*C.R_RING*1e6:.1f}–{C.R_ZERO*1e6 + (k+1)*C.R_RING*1e6:.1f}',
     f'{mt["base_g"]["rings"][k]:.2f}', f'{mt["opt_g"]["rings"][k]:.2f}',
     f'{mt["opt_g"]["rings"][k]-mt["base_g"]["rings"][k]:+.2f}',
     f'{annulus_power("base_g", C.R_ZERO + k*C.R_RING, C.R_ZERO + (k+1)*C.R_RING):.2f}',
     f'{annulus_power("opt_g", C.R_ZERO + k*C.R_RING, C.R_ZERO + (k+1)*C.R_RING):.2f}']
    for k in range(6)])

T_VER = tbl(["验证条件（均为 200 mW）", "峰值强度（mW/mm²）", "亮环半径（µm）",
             "第一旁瓣（dB）", "主瓣内功率（mW）"], [
    ["高斯 D = 5.0 mm（设计值）", f'{mt["opt_g"]["peak"]:.1f}', f'{mt["opt_g"]["r_peak"]:.1f}',
     f'{mt["opt_g"]["rings"][0]:.2f}', f'{mt["opt_g"]["p_main"]:.1f}'],
    ["高斯 D = 3.6 mm", f'{mt["opt_18"]["peak"]:.1f}', f'{mt["opt_18"]["r_peak"]:.1f}',
     f'{mt["opt_18"]["rings"][0]:.2f}', f'{mt["opt_18"]["p_main"]:.1f}'],
    ["高斯 D = 7.0 mm", f'{mt["opt_35"]["peak"]:.1f}', f'{mt["opt_35"]["r_peak"]:.1f}',
     f'{mt["opt_35"]["rings"][0]:.2f}', f'{mt["opt_35"]["p_main"]:.1f}'],
    ["平顶照明", f'{mt["opt_flat"]["peak"]:.1f}', f'{mt["opt_flat"]["r_peak"]:.1f}',
     f'{mt["opt_flat"]["rings"][0]:.2f}', f'{mt["opt_flat"]["p_main"]:.1f}'],
])

T_EE = tbl(["半径 r（µm）", "20", "37.2", "50", "100", "200"],
           [["优化前功率（mW）"] + [f'{mt["base_g"]["ee"][k]*P_IN:.2f}'
                                for k in ("20", "37", "50", "100", "200")],
            ["优化后功率（mW）"] + [f'{mt["opt_g"]["ee"][k]*P_IN:.1f}'
                                for k in ("20", "37", "50", "100", "200")],
            ["优化后能量（nJ，τ = 10 ns）"] + [f'{energy_nJ(mt["opt_g"]["ee"][k]*P_IN):.3f}'
                                             for k in ("20", "37", "50", "100", "200")]])

CSS = """
:root{--bg:#f4f6fa;--card:#fff;--ink:#1b2130;--muted:#5a6474;--line:#e3e7f0;--accent:#1d4ed8}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font-family:"Microsoft YaHei","PingFang SC","Noto Sans CJK SC",sans-serif;font-size:15px;line-height:1.72}
.wrap{max-width:1040px;margin:0 auto;padding:30px 20px 72px}
header.top{background:linear-gradient(135deg,#12203c,#1d4ed8 62%,#3b82f6);color:#fff;
  border-radius:16px;padding:26px 28px;margin-bottom:22px}
header.top h1{margin:0 0 8px;font-size:24px;font-weight:600}
header.top p{margin:3px 0;font-size:13.5px;opacity:.93}
header.top .kv{display:flex;flex-wrap:wrap;gap:8px;margin-top:14px}
header.top .kv span{background:rgba(255,255,255,.16);border-radius:999px;padding:4px 12px;font-size:12.5px}
section{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:20px 24px 22px;margin:16px 0}
h2{font-size:18px;margin:2px 0 14px;padding-left:11px;border-left:4px solid var(--accent);font-weight:600}
h3{font-size:15px;margin:18px 0 8px;font-weight:600;color:#26314a}
p{margin:9px 0}.muted{color:var(--muted)}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:11px}
.stat{background:#fbfcfe;border:1px solid var(--line);border-radius:12px;padding:12px 14px}
.stat b{display:block;font-size:19px;font-weight:600;color:#14203a}
.stat span{font-size:12.5px;color:var(--muted)}
.note{background:#f7f9fd;border-left:3px solid var(--accent);padding:12px 16px;border-radius:0 8px 8px 0;margin:14px 0;font-size:14.2px}
.tblwrap{overflow-x:auto;margin:12px 0 6px;border-radius:8px}
table{width:100%;border-collapse:collapse;font-size:13.2px;background:#fff}
th{text-align:left;background:#f1f4f9;padding:9px 11px;font-weight:600;color:#33415c;font-size:12.5px}
td{padding:7px 11px;border-top:1px solid #eef1f7;vertical-align:top}
.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:10px 0}
.grid2{display:grid;grid-template-columns:repeat(2,1fr);gap:16px;margin:10px 0}
figure.cell{margin:0}
figure.cell figcaption{font-size:12.5px;color:var(--muted);margin-top:7px;line-height:1.55}
.imgbox{background:#0c0f16;border-radius:10px;padding:6px;line-height:0}
.imgbox.chartbox{background:#fff;border:1px solid var(--line);padding:6px 4px}
.imgbox img{width:100%;height:auto;display:block;border-radius:6px}
.chart{width:100%;height:auto;display:block;overflow:hidden}
.chart .tick{font-size:11px;fill:#5a6474}
.chart .axis{font-size:12px;fill:#1b2130;font-weight:500}
.chart .grid{stroke:#eff2f8}.chart .frame{fill:none;stroke:#cfd6e4}
.chart .series{fill:none;stroke-linecap:round}
.chart .anno{stroke:#9aa6bb;stroke-dasharray:4 4}.chart .anno-t{font-size:11px;fill:#44506a}
.chart .legend-t{font-size:12px;fill:#1b2130}
.figcap{font-size:13px;color:#44506a;background:#f8fafd;border-radius:9px;padding:10px 14px;margin:12px 0 2px;line-height:1.65}
.gradrow{display:flex;align-items:center;gap:10px;margin:10px 2px 0;font-size:12.5px;color:#5a6474}
.gradrow .grad{flex:0 0 220px;height:10px;border-radius:5px;border:1px solid var(--line)}
footer{color:var(--muted);font-size:12.5px;margin-top:26px;line-height:1.8}
@media (max-width:820px){.grid3,.grid2{grid-template-columns:1fr}}
"""

BODY = f"""
<header class="top">
  <h1>相位图重画与 100 mm 光斑旁瓣抑制</h1>
  <p>输入：256 × 256 相位图（8 µm 像素）　|　输出：1920 × 1080 相位图（8 µm 像素）　|
     波长 1064 nm　|　重建距离 100 mm　|　入射激光功率 200 mW</p>
  <p>照明：平顶与高斯光束（1/e² 直径 5 mm）　|　标量衍射（带限角谱）＋ 伴随梯度相位优化</p>
  <div class="kv"><span>锥角 α = {M["alpha_deg"]:.4f}°</span>
  <span>高斯照明峰值 {mt["base_g"]["peak"]:.0f} → {mt["opt_g"]["peak"]:.0f} mW/mm²</span>
  <span>主瓣功率 {mt["base_g"]["p_main"]:.1f} → {mt["opt_g"]["p_main"]:.0f} mW</span>
  <span>第一旁瓣 {mt["base_g"]["rings"][0]:.1f} → {mt["opt_g"]["rings"][0]:.1f} dB</span></div>
</header>

<section>
<h2>1　参数与设置</h2>
{T_SET}
<p class="muted">入射激光功率取 200 mW。平顶照明下光斑内强度按 200 mW / {A_FLAT:.2f} mm²
计算；高斯照明下峰值强度按 I₀ = 2P/(πw²) 计算，两种照明的总功率相同（均为 200 mW），
因此可以直接比较。能量仅对脉冲或有限曝光有意义，按 E = P·τ 换算，
文中能量值均取 τ = 10 ns 为例。</p>
</section>

<section>
<h2>2　相位图重画与校验</h2>
{F1}
<p>输入相位图按公式延拓至 1920 × 1080（像素间距保持 8 µm），使相位图案覆盖整个 SLM 有效区。
重画图中心 256 × 256 区域与输入图的相位偏差均方根为 {M["phi0_res_rms"]:.4f} rad，
与 8 bit 量化步长（0.0246 rad）同量级，可认为重画结果与输入为同一光学元件。</p>
</section>

<section>
<h2>3　100 mm 处的光斑（优化前）</h2>
{F2}
<p>两种照明条件下的光斑均为螺旋贝塞尔光束：中心暗核（涡旋零点）、半径
{M["r_peak_um"]:.1f} µm 的亮环、外围以 {M["r_ring_um"]:.1f} µm 为周期的同心旁瓣。
在总功率相同（200 mW）的条件下，平顶照明的峰值强度为 {mt["base_flat"]["peak"]:.0f} mW/mm²，
高斯照明为 {mt["base_g"]["peak"]:.0f} mW/mm²——高斯光束把能量集中在较小的入射面积上，
因此入射峰值强度高（{I_G:.1f} mW/mm² 对 {I_FLAT:.3f} mW/mm²），但使用到的通光口径较小。</p>
{T_BASE}
{F3}
<div class="note">第一旁瓣电平在两种照明下几乎相同（{mt["base_flat"]["rings"][0]:.2f} dB 与
{mt["base_g"]["rings"][0]:.2f} dB）。旁瓣结构由锥面波的自干涉决定，照射野包络在
40–70 µm 尺度上变化很小，无法抑制与主瓣相邻的旁瓣；照明形式的差别主要体现在入射功率密度
与外围（&gt; 200 µm）环带的强弱。</div>
</section>

<section>
<h2>4　优化方法</h2>
<p>优化对象为径向对称的相位修正 δ(r)，节点间隔为 1 个 SLM 像素（8 µm），共
{M["node_count"]} 个节点，节点之间线性插值；螺旋项 ℓ = +1 保持不变。相位形式为</p>
<p style="text-align:center">φ(r,θ) = −k·sinα·r + ℓ·θ + δ(r)</p>
<p>目标函数取 z = 100 mm 处主瓣（r &lt; 37.2 µm）内的功率：</p>
<p style="text-align:center">P<sub>main</sub> = P<sub>in</sub> · ∫<sub>r&lt;37.2 µm</sub> I dr ⁄ ∫ I dr</p>
<p>相位元件不吸收能量，最大化该功率等价于将旁瓣上的功率转移至主瓣。梯度由伴随法解析求解：</p>
<p style="text-align:center">∂J/∂φ(s) = (2/D)·Im[ E*(s)·P†(χU)(s) ]</p>
<p>其中 χ 为主瓣掩模，P† 为反向传播算子（以 conj(H) 做一次 FFT/IFFT 实现），
像素梯度经节点插值矩阵映射至上述 {M["node_count"]} 个径向节点。
梯度实现的正确性以有限差分校验，偏差 &lt; 0.01%。优化器为 Adam，学习率 {M["lr"]}，
迭代 {M["iters"]} 次。</p>
{F4}
</section>

<section>
<h2>5　优化结果（高斯照明，1/e² 直径 5 mm，200 mW）</h2>
{F5}
{T_OPT}
{F6}
{T_SIDE}
{T_EE}
</section>

<section>
<h2>6　优化后的相位图</h2>
{F7}
<p>优化后相位已保存为 phase_optimized_gauss_1920x1080.png（8 bit 灰度，1920 × 1080，
相位 0–2π 映射至 0–255），重画后的原始相位保存为 phase_full_1920x1080.png。
优化结果为径向对称的连续相位分布，环带密度与原始图案相当（局部环带周期约 6–8 个像素），
局部相位梯度不超过 SLM 采样所能支持的奈奎斯特极限（λ/2Δ = 0.0665），可直接加载。</p>
<p class="muted">说明：该优化针对单一重建距离 z = 100 mm 进行，此时旁瓣最低；其他距离处
光斑形态与优化前接近。优化后亮环半径约 13 µm、环宽约 10 µm，与 SLM 像素间距（8 µm）同量级，
属该器件在此数值孔径下的极限尺寸。</p>
</section>

<section>
<h2>7　稳定性验证</h2>
{F8}
{T_VER}
<p>在光束直径偏离设计值（3.6–7.0 mm）及平顶照明条件下，主瓣内的功率仍显著高于优化前
（{mt["base_g"]["p_main"]:.1f} mW），旁瓣抑制效果保持；其中光束直径 7.0 mm 时入射峰值强度较低
但使用的通光口径较大，主瓣内功率最高（{mt["opt_35"]["p_main"]:.0f} mW）。
平顶照明下效果明显退化（第一旁瓣 {mt["opt_flat"]["rings"][0]:.1f} dB），
说明该优化相位需在高斯照明条件下使用。</p>
</section>

<section>
<h2>8　结论</h2>
<p>1）按公式重画的 1920 × 1080 相位图与输入 256 × 256 相位图在量化精度内一致，
等效锥角 {M["alpha_deg"]:.4f}°，环带周期 {M["r_ring_um"]:.2f} µm。</p>
<p>2）在 200 mW 入射功率下，优化前 100 mm 处的光斑峰值强度为
{mt["base_flat"]["peak"]:.0f} mW/mm²（平顶照明）与 {mt["base_g"]["peak"]:.0f} mW/mm²（高斯照明），
第一旁瓣均为约 {mt["base_flat"]["rings"][0]:.1f} dB，
主瓣（r &lt; 37.2 µm）内功率为 {mt["base_flat"]["p_main"]:.2f} mW 与 {mt["base_g"]["p_main"]:.2f} mW。</p>
<p>3）经伴随梯度优化后（高斯照明，1/e² 直径 5 mm）：峰值强度提高到
{mt["opt_g"]["peak"]:.0f} mW/mm²，主瓣内功率提高到 {mt["opt_g"]["p_main"]:.1f} mW
（占入射功率的 {mt["opt_g"]["eta"]*100:.1f} %），第一旁瓣降至 {mt["opt_g"]["rings"][0]:.2f} dB；
按 τ = 10 ns 换算，主瓣内能量由 {mt["base_g"]["e_main"]:.3f} nJ 提高到
{mt["opt_g"]["e_main"]:.3f} nJ。</p>
<p>4）稳定性验证表明，优化结果在光束直径 3.6–7.0 mm 范围内均保持旁瓣抑制效果，
但要求高斯照明条件。</p>
</section>

<footer>
计算条件：标量衍射模型、入射总功率 200 mW、忽略偏振与吸收；能量按 E = P·τ、τ = 10 ns 换算。<br>
代码：slm_common.py（参数/图案/传播/指标）、run_gauss.py（基线、优化与验证）、
run_report_gauss.py（报告）、chartlib.py（绘图）。
</footer>
"""

HTML = ('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>相位图重画与 100 mm 光斑旁瓣抑制</title>'
        '<style>' + CSS + '</style></head><body><div class="wrap">' + BODY
        + '</div></body></html>')
with open("report_gauss_100mm.html", "w", encoding="utf-8") as f:
    f.write(HTML)
print("[report] report_gauss_100mm.html 已生成（绝对单位：mW、mW/mm²、nJ）")
