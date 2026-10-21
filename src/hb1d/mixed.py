"""The mixed condition of Fotopoulos et al. (arXiv 2610.25104, eq. (5.4)): a necessary condition for a
heavy-ball cycle of length K in ANY dimension, hence in dimension one.

With f = m/K, M(f) = p(f) - i q(f) and A(f) = p(1 - p) - q^2,

    phi(f, Delta) = Re[(1 - M(f))^* M(f) (e^{2 pi i f Delta} - 1)] = A(f) (cos(2 pi f Delta) - 1) + q(f) sin(2 pi f Delta).

Every K-cycle (any dimension) gives weights s_m >= 0 (m = 1..K-1, proportional to |c_m|^2, not all zero) with

    sum_m s_m phi(m/K, Delta) <= 0     for Delta = 1, ..., K-1.                      (MIX_K)

So if (MIX_K) is infeasible there is no K-cycle in any dimension. By the theorem of alternatives, (MIX_K) is
infeasible iff some y >= 0 has sum_Delta y_Delta phi(m/K, Delta) > 0 for every m; such a rational y, checked with
interval arithmetic, is a rigorous certificate. Conversely (Fotopoulos et al. Prop. 5.1), feasibility of (MIX_K)
gives a cycle in dimension 2r, r = |supp s|.
"""
from fractions import Fraction as F

import numpy as np
from mpmath import iv, mp
from scipy.optimize import linprog

from .regions import _iv, _prec, sign
from .tunings import Tuning


def phi_matrix_float(t: Tuning, K: int):
    """Phi[D-1][m-1] = phi(m/K, D), floating point."""
    al, be, mu = float(t.alpha), float(t.beta), float(t.mu)
    den = (1 - mu) * al
    m = np.arange(1, K)
    f = m / K
    p = ((1 + be - mu * al) - (1 + be) * np.cos(2 * np.pi * f)) / den
    q = (1 - be) * np.sin(2 * np.pi * f) / den
    A = p * (1 - p) - q * q
    D = np.arange(1, K)[:, None]
    ang = 2 * np.pi * f[None, :] * D
    return A[None, :] * (np.cos(ang) - 1) + q[None, :] * np.sin(ang)


def phi_matrix_iv(t: Tuning, K: int, prec=256):
    with _prec(prec):
        a, b, mu = _iv(t.alpha), _iv(t.beta), _iv(t.mu)
        den = (1 - mu) * a
        out = [[None] * (K - 1) for _ in range(K - 1)]
        for m in range(1, K):
            ang = 2 * iv.pi * m / K
            p = ((1 + b - mu * a) - (1 + b) * iv.cos(ang)) / den
            q = (1 - b) * iv.sin(ang) / den
            A = p * (1 - p) - q * q
            for D in range(1, K):
                out[D - 1][m - 1] = A * (iv.cos(ang * D) - 1) + q * iv.sin(ang * D)
        return out


def _rationalize(v, den=10 ** 12):
    return [max(F(0), F(float(x)).limit_denominator(den)) for x in v]


def mixed(t: Tuning, K: int, prec=256):
    """Decide (MIX_K) for a rational tuning.

    Returns dict(feasible=True/False/None, margin=float, cert=...). feasible=False comes with a rational y whose
    certified products sum_D y_D phi(m/K, D) are > 0 for every m (rigorous). feasible=True comes with rational
    weights s whose certified constraint values are all < 0 (strict feasibility, rigorous); boundary cases with
    zero margin return feasible=None."""
    Phi = phi_matrix_float(t, K)
    n = K - 1
    # alternative: max tau s.t. Phi^T y >= tau, sum y = 1, y >= 0
    c = np.zeros(n + 1)
    c[-1] = -1
    A_ub = np.hstack([-Phi.T, np.ones((n, 1))])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(n), A_eq=[[1] * n + [0]], b_eq=[1],
                  bounds=[(0, None)] * n + [(None, None)], method="highs")
    tau = -res.fun
    if tau > 1e-12:
        y = _rationalize(res.x[:n])
        P = phi_matrix_iv(t, K, prec)
        with _prec(prec):
            vals = [sum((_iv(y[D]) * P[D][m] for D in range(n) if y[D]), iv.mpf(0)) for m in range(n)]
        if all(sign(v) == 1 for v in vals):
            return dict(feasible=False, margin=tau, cert=y)
        return dict(feasible=None, margin=tau, note="infeasibility certificate failed rigorous check")
    # primal: max sigma s.t. Phi s <= -sigma, sum s = 1, s >= 0
    A_ub = np.hstack([Phi, np.ones((n, 1))])
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(n), A_eq=[[1] * n + [0]], b_eq=[1],
                  bounds=[(0, None)] * n + [(None, None)], method="highs")
    sigma = -res.fun
    if sigma > 1e-12:
        s = _rationalize(res.x[:n])
        P = phi_matrix_iv(t, K, prec)
        with _prec(prec):
            vals = [sum((_iv(s[m]) * P[D][m] for m in range(n) if s[m]), iv.mpf(0)) for D in range(n)]
        if all(sign(v) == -1 for v in vals):
            support = [m + 1 for m in range(n) if s[m] > 0]
            return dict(feasible=True, margin=sigma, cert=s, support=support)
        return dict(feasible=None, margin=sigma, note="feasibility certificate failed rigorous check")
    return dict(feasible=None, margin=max(tau, sigma), note="boundary: zero margin on both sides")


def check_weights(t: Tuning, K: int, weights: dict, prec=256):
    """Rigorously check given rational weights {m: s_m} strictly satisfy (MIX_K). Returns certified max value."""
    P = phi_matrix_iv(t, K, prec)
    with _prec(prec):
        vals = [sum((_iv(w) * P[D][m - 1] for m, w in weights.items()), iv.mpf(0)) for D in range(K - 1)]
    return all(sign(v) == -1 for v in vals), max(float(mp.mpf(v.b)) for v in vals)


def fan_program_slacks(t: Tuning, k: int, nu: dict):
    """Pasechnyuk-Vilensky Thm 'Fan tongues'(a) at exponent t = 0 for a k-periodic spectrum nu = {j: weight}
    (frequencies omega_j = 2 pi j/k): with w_j = |g - 2is|(omega_j) nu_j, where |g - 2is| = |a(z) b(z)|/B^2 on
    |z| = 1, returns the slacks [sum_j w_j Lambda_n(omega_j) for n = 1..k-1] + [sum_j w_j cos theta_j]; the program
    is satisfied iff all are >= 0. Lambda_n(omega) = cos Theta - cos(n omega - Theta). Floating point."""
    import cmath
    import math

    from .regions import theta_def
    al, be, mu = float(t.alpha), float(t.beta), float(t.mu)
    B = al * (1 - mu) / 2
    ws, ths, oms = [], [], []
    for j, v in nu.items():
        om = 2 * math.pi * j / k
        z = cmath.exp(1j * om)
        a = z * z - (1 + be - al * mu) * z + be
        b = z * z - (1 + be - al) * z + be
        ws.append(float(v) * abs(a) * abs(b) / B ** 2)
        ths.append(theta_def(t, om))
        oms.append(om)
    lags = [sum(w * (math.cos(th) - math.cos(n * om - th)) for w, th, om in zip(ws, ths, oms)) for n in range(1, k)]
    return lags + [sum(w * math.cos(th) for w, th in zip(ws, ths))]
