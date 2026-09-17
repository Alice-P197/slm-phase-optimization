# -*- coding: utf-8 -*-
"""
径向相位优化模块：目标函数、伴随梯度、Adam 迭代。

优化变量
--------
径向相位修正 delta(r)，共 1103 个节点（间隔 1 个 SLM 像素 = 8 µm），线性插值；
螺旋项 ell*theta 固定不变，相位 phi(r,theta) = -k*sin(alpha)*r + ell*theta + delta(r)。

目标函数
--------
z = 100 mm 处主瓣（r < 37.18 µm）内的能量占比 J = sum(CHI*I)/sum(I)。
相位元件不吸收能量，最大化 J 等价于把旁瓣上的功率转移回主瓣。

梯度
----
伴随法：dJ/dphi(s) = (2/D) * Im[ conj(E(s)) * P_dagger(CHI*U)(s) ]。
P_dagger 即反向传播（传递函数取共轭做一次 FFT/IFFT），与神经网络反向传播同构。
"""

import time
import numpy as np
import slm_common as C

PAD = 64
Z = 0.100
R_C = C.R_ZERO                   # 主瓣半径 37.18 µm
W_MAIN = 2.5e-3                  # 1/e² 强度半径，对应光束直径 5 mm
NODE_STEP = 8e-6
ITERS = 200
LR = 0.05
B1, B2, EPS = 0.9, 0.999, 1e-8   # Adam 超参数


class RadialPhaseOptimizer:
    """螺旋锥透镜相位图的径向相位优化器。"""

    def __init__(self, gauss_w=W_MAIN, z=Z, r_core=R_C, node_step=NODE_STEP,
                 pad=PAD, seed=0):
        self.gauss_w = gauss_w
        self.z = z
        self.r_core = r_core
        self.node_step = node_step
        self.pad = pad
        self.seed = seed

        # ---- 由输入 256x256 相位图标定参数，并按公式重画整幅 SLM ----
        self.phi0, self.res_rms, self.phi_in = C.fit_phi0()
        self.X, self.Y = C.grid_coords()
        self.R = np.hypot(self.X, self.Y)
        self.TH = np.arctan2(self.Y, self.X)
        self.BASE = -C.K * C.SIN_A * self.R
        self.phase_full = np.mod(self.BASE + C.ELL * self.TH + self.phi0, C.FULL)

        # ---- 径向节点参数化（线性插值） ----
        nodes = np.arange(0.0, C.R_CORNER + self.node_step, self.node_step)
        self.NN = nodes.size
        tt = np.clip(self.R / self.node_step, 0, self.NN - 1.001)
        self.ni = np.floor(tt).astype(int)
        self.fr = tt - self.ni

        # ---- 补零传播网格与主瓣掩模 ----
        HH, WW = C.SLM_H + 2 * self.pad, C.SLM_W + 2 * self.pad
        self.HH, self.WW = HH, WW
        self.H = C.transfer_function(self.z, WW, HH)
        self.HI = np.conj(self.H)
        yy, xx = np.mgrid[0:HH, 0:WW].astype(float)
        self.CEN = (C.SLM_W / 2.0 + self.pad - 0.5, C.SLM_H / 2.0 + self.pad - 0.5)
        rr = np.hypot((xx - self.CEN[0]) * C.DX, (yy - self.CEN[1]) * C.DX)
        self.CHI = (rr <= self.r_core).astype(float)

    # -------------------------------------------------------------- 场与传播
    def field(self, delta, gauss_w):
        phi = (self.BASE + C.ELL * self.TH
               + delta[self.ni] * (1 - self.fr) + delta[self.ni + 1] * self.fr)
        return C.incident_amplitude(self.X, self.Y, gauss_w) * np.exp(1j * phi)

    def fwd(self, E):
        Ep = np.zeros((self.HH, self.WW), np.complex128)
        Ep[self.pad:self.pad + C.SLM_H, self.pad:self.pad + C.SLM_W] = E
        return np.fft.ifft2(np.fft.fft2(Ep) * self.H)

    def bwd(self, V):
        U = np.fft.ifft2(np.fft.fft2(V) * self.HI)
        return U[self.pad:self.pad + C.SLM_H, self.pad:self.pad + C.SLM_W]

    # ------------------------------------------------------ 目标函数 + 梯度
    def objective_and_gradient(self, delta, want_grad=True):
        E = self.field(delta, self.gauss_w)
        U = self.fwd(E)
        I = np.abs(U) ** 2
        D = float(I.sum())
        J = float((self.CHI * I).sum() / D)
        if not want_grad:
            return J, None
        g = self.bwd(self.CHI * U / D)
        gp = 2.0 * np.imag(np.conj(E) * g)                 # 像素梯度 dJ/dphi
        gn = np.bincount(self.ni.ravel(), weights=(gp * (1 - self.fr)).ravel(),
                         minlength=self.NN)
        gn += np.bincount((self.ni + 1).ravel(), weights=(gp * self.fr).ravel(),
                          minlength=self.NN)               # 映射到径向节点
        return J, gn

    # ---------------------------------------------------------- 梯度有限差分校验
    def gradient_check(self, delta):
        J0, g0 = self.objective_and_gradient(delta)
        rng = np.random.default_rng(self.seed)
        v = rng.normal(size=self.NN)
        v /= np.linalg.norm(v)
        h = 2e-4
        fd = (self.objective_and_gradient(delta + h * v, False)[0]
              - self.objective_and_gradient(delta - h * v, False)[0]) / (2 * h)
        an = float(np.dot(g0, v))
        return dict(J0=float(J0), finite_diff=float(fd), adjoint=an,
                    rel_err_pct=abs(fd - an) / abs(fd) * 100)

    # ---------------------------------------------------------------- Adam 迭代
    def run(self, delta0=None, iters=ITERS, lr=LR, verbose=True):
        delta = np.zeros(self.NN) if delta0 is None else np.asarray(delta0).copy()
        chk = self.gradient_check(delta)
        if verbose:
            print(f"    [optimize] 初始主瓣能量 {chk['J0']*100:.3f}%；"
                  f"伴随梯度与有限差分偏差 {chk['rel_err_pct']:.3f}%", flush=True)
        m_ad = np.zeros(self.NN)
        v_ad = np.zeros(self.NN)
        hist = []
        best_J, best_d = chk["J0"], delta.copy()
        t0 = time.time()
        for it in range(1, iters + 1):
            J, g = self.objective_and_gradient(delta)
            m_ad = B1 * m_ad + (1 - B1) * g
            v_ad = B2 * v_ad + (1 - B2) * g ** 2
            delta = (delta + lr * (m_ad / (1 - B1 ** it))
                     / (np.sqrt(v_ad / (1 - B2 ** it)) + EPS))
            hist.append(float(J))
            if J > best_J:
                best_J, best_d = float(J), delta.copy()
            if verbose and (it % 20 == 0 or it == 1):
                print(f"    [optimize] it {it:4d}  主瓣能量 {J*100:6.3f}%  "
                      f"t={time.time()-t0:.0f}s", flush=True)
        return best_d, np.asarray(hist), best_J, chk
