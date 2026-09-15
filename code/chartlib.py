# -*- coding: utf-8 -*-
"""Small chart/figure helpers (pure numpy + PIL + inline SVG)."""

import base64
import io
import numpy as np
from PIL import Image

PAL = ["#1d4ed8", "#d97706", "#15803d", "#7c3aed", "#dc2626", "#0891b2"]
INK, MUTED, GRID = "#1b2130", "#5a6474", "#dfe4ee"

_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def sup10(k):
    return "10" + str(int(k)).translate(_SUP)


def _hex(c):
    return "#%02x%02x%02x" % tuple(int(round(v * 255)) for v in c)


MAGMA = [(0.00, (0.001, 0.000, 0.014)), (0.13, (0.110, 0.066, 0.270)),
         (0.25, (0.317, 0.071, 0.485)), (0.38, (0.510, 0.145, 0.510)),
         (0.50, (0.712, 0.212, 0.475)), (0.63, (0.897, 0.316, 0.392)),
         (0.75, (0.984, 0.529, 0.382)), (0.88, (0.996, 0.761, 0.529)),
         (1.00, (0.988, 0.992, 0.749))]
DIVERGE = [(0.0, (0.129, 0.400, 0.674)), (0.25, (0.569, 0.749, 0.898)),
           (0.5, (0.968, 0.968, 0.968)), (0.75, (0.937, 0.545, 0.396)),
           (1.0, (0.698, 0.094, 0.169))]
GRAY = [(0.0, (0.08, 0.08, 0.09)), (1.0, (0.98, 0.98, 0.98))]


def lut(stops, n=256):
    stops = sorted(stops)
    p = np.array([s[0] for s in stops])
    c = np.array([s[1] for s in stops])
    t = np.linspace(0, 1, n)
    return np.stack([np.interp(t, p, c[:, i]) for i in range(3)], 1)


def hsv_lut(n=256):
    t = np.linspace(0, 1, n, endpoint=False)
    h = t * 6
    i = np.floor(h).astype(int) % 6
    q = 1 - (h - i)
    one, zero = np.ones(n), np.zeros(n)
    sectors = [(one, q, zero), (q, one, zero), (zero, one, q),
               (zero, q, one), (q, zero, one), (one, zero, q)]
    r, g, b = np.zeros(n), np.zeros(n), np.zeros(n)
    for k, (rv, gv, bv) in enumerate(sectors):
        m = i == k
        r[m], g[m], b[m] = rv[m], gv[m], bv[m]
    return np.stack([r, g, b], 1)


LUT_MAG = lut(MAGMA)
LUT_DIV = lut(DIVERGE)
LUT_HSV = hsv_lut()
LUT_GRAY = lut(GRAY)


def norm(Z, log=False, vmin=None, vmax=None, gamma=1.0):
    Z = np.asarray(Z, float)
    if log:
        Z = np.log10(np.maximum(Z, 1e-12))
    lo = np.nanmin(Z) if vmin is None else vmin
    hi = np.nanmax(Z) if vmax is None else vmax
    return np.clip((Z - lo) / (hi - lo + 1e-30), 0, 1) ** gamma


def png_uri(arr01, cmap, size=None):
    idx = np.clip((arr01 * 255).round().astype(int), 0, 255)
    rgb = (cmap[idx] * 255).astype(np.uint8)
    im = Image.fromarray(rgb, "RGB")
    if size:
        im = im.resize(size, Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def crop(A, half_px, centre=None):
    c = A.shape[0] // 2 if centre is None else int(round(centre))
    return A[c - half_px:c + half_px, c - half_px:c + half_px]


def css_gradient(stops, direction="to right"):
    # 色标必须写成 "<颜色> <位置>"（颜色在前），否则 Chromium 会判定整个声明无效
    p = [f"{_hex(s[1])} {s[0]*100:.1f}%" for s in sorted(stops)]
    return f"linear-gradient({direction}, " + ", ".join(p) + ")"


def nice_ticks(lo, hi, n=6):
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        hi = lo + 1
    raw = (hi - lo) / max(n, 2)
    mag = 10 ** np.floor(np.log10(raw))
    step = mag * 10
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            step = m * mag
            break
    return np.arange(np.ceil(lo / step) * step, hi + step * 0.4, step), step


def fmt(v, step):
    if step >= 1:
        return f"{v:.0f}"
    if step >= 0.1:
        return f"{v:.1f}"
    if step >= 0.01:
        return f"{v:.2f}"
    return f"{v:g}"


def chart(series, xlabel, ylabel, w=700, h=330, logy=False, xlim=None, ylim=None,
          annotations=(), nxt=6, nyt=6, legpos="ne"):
    ml, mr, mt, mb = 78, 16, 16, 48
    px0, px1, py0, py1 = ml, w - mr, mt, h - mb
    xs = np.concatenate([np.asarray(s["x"], float) for s in series])
    ys = np.concatenate([np.asarray(s["y"], float) for s in series])
    if xlim is None:
        xlim = (float(np.nanmin(xs)), float(np.nanmax(xs)))
    if logy:
        yy = np.log10(np.maximum(ys, 1e-12))
        ylim_ = (float(np.nanmin(yy)), float(np.nanmax(yy))) if ylim is None else ylim
    elif ylim is None:
        lo, hi = float(np.nanmin(ys)), float(np.nanmax(ys))
        if lo > 0:
            lo = 0.0
        if hi < 0:
            hi = 0.0
        pad = (hi - lo) * 0.06
        ylim_ = (lo - pad if lo < 0 else 0.0, hi + pad)
    else:
        ylim_ = ylim

    def sx(v):
        return px0 + (np.asarray(v, float) - xlim[0]) / (xlim[1] - xlim[0]) * (px1 - px0)

    def sy(v):
        v = np.log10(np.maximum(np.asarray(v, float), 1e-12)) if logy else np.asarray(v, float)
        return py1 - (v - ylim_[0]) / (ylim_[1] - ylim_[0]) * (py1 - py0)

    xt, xstep = nice_ticks(xlim[0], xlim[1], nxt)
    if logy:
        lo, hi = ylim_
        if hi - lo <= 2.2:
            cand = [k * 10.0 ** d for d in range(int(np.floor(lo)) - 1, int(np.ceil(hi)) + 1)
                    for k in (1, 2, 5)]
            yt = np.array([np.log10(c) for c in cand if lo - 1e-9 <= np.log10(c) <= hi + 1e-9])
            ylab = [f"{10**v:g}" for v in yt]
        else:
            yt = np.arange(np.floor(lo), np.ceil(hi) + 0.001)
            yt = yt[(yt >= lo - 1e-9) & (yt <= hi + 1e-9)]
            ylab = [sup10(v) for v in yt]
    else:
        yt, ystep = nice_ticks(ylim_[0], ylim_[1], nyt)
        ylab = [fmt(v, ystep) for v in yt]

    o = [f'<svg class="chart" viewBox="0 0 {w} {h}" role="img" '
         f'aria-label="{xlabel} / {ylabel}" preserveAspectRatio="xMidYMid meet">']
    for v in yt:
        y = float(sy(10 ** v if logy else v))
        o.append(f'<line class="grid" x1="{px0}" x2="{px1}" y1="{y:.1f}" y2="{y:.1f}"/>')
    o.append(f'<rect class="frame" x="{px0}" y="{py0}" width="{px1-px0}" height="{py1-py0}"/>')
    for v, lb in zip(yt, ylab):
        y = float(sy(10 ** v if logy else v))
        o.append(f'<text class="tick" x="{px0-8}" y="{y+4:.1f}" text-anchor="end">{lb}</text>')
    for v in xt:
        if v < xlim[0] - 1e-12 or v > xlim[1] + 1e-12:
            continue
        x = float(sx(v))
        o.append(f'<line class="tick" x1="{x:.1f}" x2="{x:.1f}" y1="{py1}" y2="{py1+5}"/>')
        o.append(f'<text class="tick" x="{x:.1f}" y="{py1+20}" text-anchor="middle">'
                 f'{fmt(v, xstep)}</text>')
    o.append(f'<text class="axis" x="{(px0+px1)/2:.0f}" y="{h-8}" text-anchor="middle">{xlabel}</text>')
    o.append(f'<text class="axis" x="18" y="{(py0+py1)/2:.0f}" text-anchor="middle" '
             f'transform="rotate(-90 18 {(py0+py1)/2:.0f})">{ylabel}</text>')
    for s in series:
        X, Y = np.asarray(s["x"], float), np.asarray(s["y"], float)
        m = np.isfinite(X) & np.isfinite(Y)
        if logy:
            m &= Y > 0
        X, Y = X[m], Y[m]
        if X.size == 0:
            continue
        if X[0] > X[-1]:
            o2 = np.argsort(X)
            X, Y = X[o2], Y[o2]
        px, py = sx(X), sy(Y)
        d = "M " + " L ".join(f"{a:.1f},{b:.1f}" for a, b in zip(px, py))
        style = f'stroke="{s.get("color", PAL[0])}"'
        if s.get("dash"):
            style += f' stroke-dasharray="{s["dash"]}"'
        o.append(f'<path class="series" d="{d}" {style} stroke-width="{s.get("width", 2)}"/>')
    for a in annotations:
        if a[0] == "vline":
            x = float(sx(a[1]))
            if not (px0 - 1 <= x <= px1 + 1):
                continue
            dy = float(a[3]) if len(a) > 3 else 0.0
            anchor = "end" if x > (px0 + px1) / 2 else "start"
            tx = x - 4 if anchor == "end" else x + 4
            o.append(f'<line class="anno" x1="{x:.1f}" x2="{x:.1f}" y1="{py0}" y2="{py1}"/>')
            o.append(f'<text class="anno-t" x="{tx:.1f}" y="{py0+12+dy:.1f}" '
                     f'text-anchor="{anchor}">{a[2]}</text>')
        elif a[0] == "hline":
            y = float(sy(a[1]))
            if not (py0 - 1 <= y <= py1 + 1):
                continue
            o.append(f'<line class="anno" x1="{px0}" x2="{px1}" y1="{y:.1f}" y2="{y:.1f}"/>')
            o.append(f'<text class="anno-t" x="{px0+6}" y="{y-5:.1f}">{a[2]}</text>')
        elif a[0] == "dot":
            x, y = float(sx(a[1])), float(sy(a[2]))
            dy = float(a[4]) if len(a) > 4 else 0.0
            o.append(f'<circle class="dot" cx="{x:.1f}" cy="{y:.1f}" r="3.5"/>')
            o.append(f'<text class="anno-t" x="{x+7:.1f}" y="{y-7+dy:.1f}">{a[3]}</text>')
    lab = [s for s in series if s.get("label")]
    if lab:
        if legpos == "ne":
            lx, ly = px1 - 210, py0 + 14
        elif legpos == "sw":
            lx, ly = px0 + 10, py1 - 20 * len(lab) - 6
        else:
            lx, ly = px0 + 8, py0 + 14
        o.append('<g class="legend">')
        for i, s in enumerate(lab):
            yy = ly + i * 17
            o.append(f'<line x1="{lx}" x2="{lx+22}" y1="{yy-4}" y2="{yy-4}" '
                     f'stroke="{s.get("color", PAL[0])}" stroke-width="2.5"'
                     + (f' stroke-dasharray="{s["dash"]}"' if s.get("dash") else "") + '/>')
            o.append(f'<text class="legend-t" x="{lx+28}" y="{yy}">{s["label"]}</text>')
        o.append('</g>')
    o.append("</svg>")
    return "".join(o)


def svg_img(uri, xlab, ylab, xr, yr, w=720, h=340, xt=None, yt=None):
    ml, mr, mt, mb = 80, 16, 16, 48
    px0, px1, py0, py1 = ml, w - mr, mt, h - mb
    o = [f'<svg class="chart" viewBox="0 0 {w} {h}" role="img" aria-label="{xlab} / {ylab}">']
    o.append(f'<image href="{uri}" x="{px0}" y="{py0}" width="{px1-px0}" height="{py1-py0}" '
             f'preserveAspectRatio="none"/>')
    o.append(f'<rect class="frame" x="{px0}" y="{py0}" width="{px1-px0}" height="{py1-py0}"/>')
    if xt is None:
        xt, _ = nice_ticks(*xr, 6)
    for v in xt:
        x = px0 + (v - xr[0]) / (xr[1] - xr[0]) * (px1 - px0)
        o.append(f'<line class="tick" x1="{x:.1f}" x2="{x:.1f}" y1="{py1}" y2="{py1+5}"/>')
        o.append(f'<text class="tick" x="{x:.1f}" y="{py1+20}" text-anchor="middle">{v:g}</text>')
    for v in yt:
        y = py1 - (v - yr[0]) / (yr[1] - yr[0]) * (py1 - py0)
        o.append(f'<line class="tick" x1="{px0-5}" x2="{px0}" y1="{y:.1f}" y2="{y:.1f}"/>')
        o.append(f'<text class="tick" x="{px0-8}" y="{y+4:.1f}" text-anchor="end">{v:g}</text>')
    o.append(f'<text class="axis" x="{(px0+px1)/2:.0f}" y="{h-8}" text-anchor="middle">{xlab}</text>')
    o.append(f'<text class="axis" x="18" y="{(py0+py1)/2:.0f}" text-anchor="middle" '
             f'transform="rotate(-90 18 {(py0+py1)/2:.0f})">{ylab}</text>')
    o.append("</svg>")
    return "".join(o)


def tbl(head, rows):
    h = "".join(f"<th>{x}</th>" for x in head)
    b = "".join("<tr>" + "".join(f"<td>{x}</td>" for x in r) + "</tr>" for r in rows)
    return (f'<div class="tblwrap"><table><thead><tr>{h}</tr></thead>'
            f'<tbody>{b}</tbody></table></div>')


def figcap(text):
    return f'<p class="figcap">{text}</p>'
