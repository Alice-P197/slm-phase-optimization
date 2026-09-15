# -*- coding: utf-8 -*-
"""
公共模块：物理参数、相位图案生成、衍射传播、指标计算。

物理设定
--------
SLM            : 1920 x 1080 像素，像素间距 8 µm  ->  15.36 mm x 8.64 mm
输入相位图     : 256 x 256，像素间距 8 µm（等效锥角由此标定）
波长           : 1064 nm
相位形式       : phi(r,theta) = -k*sin(alpha)*r + ell*theta + phi0   (mod 2π)
                 即"螺旋锥透镜":锥面相位 + 拓扑荷 ell 的螺旋相位
"""

import numpy as np
from PIL import Image

# ------------------------------------------------------------------ 常量
LAM = 1.064e-6                 # 波长 [m]
DX = 8e-6                      # SLM 像素间距 [m]
SLM_W, SLM_H = 1920, 1080      # SLM 像素数
AP_W, AP_H = SLM_W * DX, SLM_H * DX        # 15.36 mm x 8.64 mm

SLOPE = 0.824481               # 输入相位图测得的径向相位斜率 [rad/pixel @ 8µm]
ELL = 1.0                      # 拓扑荷
SIN_A = SLOPE * LAM / (2 * np.pi * DX)     # = 0.0174529
ALPHA = float(np.arcsin(SIN_A))            # = 1.0000 deg
K = 2 * np.pi / LAM

R_HALF_V = AP_H / 2            # 4.32 mm   (垂直半宽)
R_HALF_H = AP_W / 2            # 7.68 mm   (水平半宽)
R_CORNER = float(np.hypot(R_HALF_H, R_HALF_V))   # 8.83 mm
Z_MAX = R_HALF_V / np.tan(ALPHA)                 # 247.5 mm (按垂直半宽)

# 理想贝塞尔-涡旋光斑的特征尺寸
R_PEAK = 1.8412 * LAM / (2 * np.pi * SIN_A)      # 17.87 µm
R_ZERO = 3.8317 * LAM / (2 * np.pi * SIN_A)      # 37.18 µm
R_RING = np.pi * LAM / (2 * np.pi * SIN_A)       # 30.49 µm
FULL = 2 * np.pi


# ------------------------------------------------------------- 图案生成
def grid_coords(w=SLM_W, h=SLM_H, dx=DX):
    """SLM 平面的坐标网格（原点在像素阵列中心）。"""
    x = (np.arange(w) - (w - 1) / 2.0) * dx
    y = (np.arange(h) - (h - 1) / 2.0) * dx
    X, Y = np.meshgrid(x, y)                     # Y: h x w
    return X, Y


def spiral_axicon_phase(x=None, y=None, sin_a=SIN_A, ell=ELL, phi0=0.0,
                        w=SLM_W, h=SLM_H, dx=DX):
    """按公式生成螺旋锥透镜相位，铺满整幅 SLM。返回 (相位, 半径, 方位角)。"""
    if x is None:
        X, Y = grid_coords(w, h, dx)
    else:
        X, Y = x, y
    R = np.hypot(X, Y)
    TH = np.arctan2(Y, X)
    phi = -K * sin_a * R + ell * TH + phi0
    return np.mod(phi, FULL), R, TH


def fit_phi0(path=r"C:\Users\Fischer\Desktop\mask_phase_hologram.png"):
    """用输入的 256x256 相位图标定常数相位偏移 phi0，并给出残差。"""
    a = np.array(Image.open(path))
    if a.ndim == 3:
        a = a[:, :, 0]
    phi_in = a.astype(float) / 255.0 * FULL
    n = phi_in.shape[0]
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    R = np.hypot(xx - c, yy - c) * DX
    TH = np.arctan2(yy - c, xx - c)
    model = -K * SIN_A * R + ELL * TH
    d = np.mod(phi_in - model, FULL)
    m = R > 6 * DX
    phi0 = float(np.angle(np.mean(np.exp(1j * d[m]))))            # 圆均值
    res = np.mod(d - phi0 + np.pi, FULL) - np.pi
    return phi0, float(np.sqrt((res[m] ** 2).mean())), phi_in


# ------------------------------------------------------------- 传播
def transfer_function(z, w, h, dx=DX, lam=LAM, band_limit=True):
    """带限角谱传递函数: H = exp(i*2πz/λ*sqrt(1-(λfx)²-(λfy)²))。"""
    fx = np.fft.fftfreq(w, dx)
    fy = np.fft.fftfreq(h, dx)
    FX, FY = np.meshgrid(fx, fy)
    arg = 1.0 - (lam * FX) ** 2 - (lam * FY) ** 2
    H = np.zeros_like(arg, dtype=np.complex128)
    m = arg > 0
    H[m] = np.exp(2j * np.pi * z / lam * np.sqrt(arg[m]))
    if band_limit:
        flim_x = 1.0 / (lam * np.sqrt((2 * (1.0 / (w * dx)) * z) ** 2 + 1.0))
        flim_y = 1.0 / (lam * np.sqrt((2 * (1.0 / (h * dx)) * z) ** 2 + 1.0))
        H[np.abs(FX) > flim_x] = 0
        H[np.abs(FY) > flim_y] = 0
    return H


def pad_field(E, w, h, pad=64):
    """把 SLM 上的场嵌入更大的零场（避免 FFT 环绕）。"""
    out = np.zeros((h + 2 * pad, w + 2 * pad), np.complex128)
    out[pad:pad + h, pad:pad + w] = E
    return out


def propagate(E, z, dx=DX):
    """标量衍射传播（带限角谱）: 输入/输出均为同一网格上的复振幅。"""
    h, w = E.shape
    H = transfer_function(z, w, h, dx)
    return np.fft.ifft2(np.fft.fft2(E) * H)


def propagate_windows(Eslm, z, dx=DX, pad=64):
    """把 SLM 场补零后传播，返回补零网格上的强度与中心坐标。"""
    h, w = Eslm.shape
    E = pad_field(Eslm, w, h, pad)
    U = propagate(E, z, dx)
    I = np.abs(U) ** 2
    return I, (w / 2.0 + pad - 0.5, h / 2.0 + pad - 0.5)


# ------------------------------------------------------------- 指标
def radius_map(I, centre, dx=DX):
    h, w = I.shape
    y, x = np.mgrid[0:h, 0:w]
    return np.hypot((x - centre[0]) * dx, (y - centre[1]) * dx)


def radial_profile(I, centre, dx=DX, rmax=None, nbins=1200):
    """环带平均径向剖面；空环带用相邻有效值线性插值填补。"""
    r = radius_map(I, centre, dx)
    if rmax is None:
        rmax = float(r.max())
    edges = np.linspace(0, rmax, nbins + 1)
    idx = np.digitize(r.ravel(), edges) - 1
    ok = (idx >= 0) & (idx < nbins)
    s = np.bincount(idx[ok], weights=I.ravel()[ok], minlength=nbins)
    n = np.bincount(idx[ok], minlength=nbins)
    rc = 0.5 * (edges[1:] + edges[:-1])
    p = np.where(n > 0, s / np.maximum(n, 1), np.nan)
    good = ~np.isnan(p)
    if good.sum() < 3:
        return rc, np.zeros_like(rc)
    return rc, np.interp(rc, rc[good], p[good])


def annulus_max(I, centre, r1, r2, dx=DX):
    r = radius_map(I, centre, dx)
    m = (r >= r1) & (r < r2)
    if not np.any(m):
        return 0.0
    return float(I[m].max())


def encircled(I, centre, radii, dx=DX):
    r = radius_map(I, centre, dx)
    tot = float(I.sum())
    return {q: float(I[r <= q].sum() / tot) for q in radii}


def sidelobe_table(I, centre, rc=R_ZERO, spacing=R_RING, nring=6, dx=DX):
    """主瓣峰值 + 各旁瓣环带的峰值电平（dB，相对主瓣峰值）。"""
    r = radius_map(I, centre, dx)
    main = float(I[r <= rc].max())
    main_r = float(r[r <= rc][np.argmax(I[r <= rc])])
    rows = []
    for k in range(nring):
        lo, hi = rc + k * spacing, rc + (k + 1) * spacing
        v = annulus_max(I, centre, lo, hi, dx)
        rows.append((k + 1, lo, hi, v, 10 * np.log10(v / main) if v > 0 else -np.inf))
    return main, main_r, rows


def core_metrics(I, centre, dx=DX, rc=R_ZERO):
    """主瓣峰值、峰值半径、第一零点半径、主瓣/旁瓣能量占比。"""
    r = radius_map(I, centre, dx)
    m = r <= rc
    main = float(I[m].max())
    main_r = float(r[m][np.argmax(I[m])])
    tot = float(I.sum())
    eta_in = float(I[m].sum() / tot)
    rc_p, prof = radial_profile(I, centre, dx, rmax=rc * 3, nbins=900)
    ipk = int(np.argmax(prof))
    k = ipk
    while k < prof.size - 1 and prof[k + 1] < prof[k]:
        k += 1
    return dict(peak=main, r_peak=main_r, r_zero=float(rc_p[k]),
                eta_main=eta_in, eta_side=1.0 - eta_in)


def phase_to_png(phase, path, bits=8):
    """把 0..2π 的相位保存为 SLM 用的灰度 PNG。"""
    v = np.mod(phase, FULL) / FULL * (2 ** bits - 1)
    dt = np.uint8 if bits == 8 else np.uint16
    Image.fromarray(np.round(v).clip(0, 2 ** bits - 1).astype(dt)).save(path)


# ------------------------------------------------------------- 照明
def incident_amplitude(X, Y, gauss_w=None):
    """入射光振幅分布。

    gauss_w = None : 单位振幅平面波（平顶照明）
    gauss_w = w    : 高斯光束，振幅 A(r) = exp(-r²/w²)，
                     即强度半宽（1/e² 半径）为 w，峰值归一化为 1。
    """
    if gauss_w is None:
        return np.ones_like(X)
    return np.exp(-(X ** 2 + Y ** 2) / gauss_w ** 2)


def transmission(gauss_w=None, w=SLM_W, h=SLM_H, dx=DX):
    """入射总功率（相对峰值强度），用于效率归一化。"""
    X, Y = grid_coords(w, h, dx)
    A = incident_amplitude(X, Y, gauss_w)
    return float((A ** 2).sum()),
