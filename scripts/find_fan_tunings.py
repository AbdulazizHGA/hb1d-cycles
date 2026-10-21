"""Find rational tunings at kappa = 100 strictly inside the first-level fan tongues 2/(2K+1) of Pasechnyuk-Vilensky,
with a certified distance from Omega_aug.

For a gap I_K the 2/(2K+1) fan tongue sits at the top of the gap, next to the (K+1)-tongue. We scan the gap in
unrolled coordinates (rho, delta/delta_top), then maximize min(mixed margin at level k, -max pure-mode value) with
Nelder-Mead, round to small denominators and certify both properties rigorously (hb1d.mixed.mixed,
hb1d.regions.omega_aug_report). Writes results/fan_tunings.json.

    python scripts/find_fan_tunings.py
"""
import json
import math
import os
import sys
from fractions import Fraction as F

import numpy as np
from scipy.optimize import brentq, linprog, minimize

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from hb1d.geometry import lag_one, tuning_to_unrolled, unrolled_to_tuning  # noqa: E402
from hb1d.mixed import mixed, phi_matrix_float  # noqa: E402
from hb1d.regions import omega_aug_report  # noqa: E402
from hb1d.tunings import Tuning  # noqa: E402

KAPPA = F(100)
MU = 1 / float(KAPPA)


def tun(a, b):
    return Tuning("candidate", F(a), F(b), KAPPA)


def mixed_margin(a, b, k):
    Phi = phi_matrix_float(tun(a, b), k)
    n = k - 1
    res = linprog(np.r_[np.zeros(n), -1], A_ub=np.hstack([Phi, np.ones((n, 1))]), b_ub=np.zeros(n),
                  A_eq=[[1] * n + [0]], b_eq=[1], bounds=[(0, None)] * n + [(None, None)])
    return -res.fun


def max_pure_mode(a, b, Kmax=40):
    """max over primitive modes K <= Kmax of A - q cot(pi/K) (float); negative = outside these Omega_{K,m}."""
    best = -math.inf
    den = (1 - MU) * a
    for K in range(3, Kmax + 1):
        for m in range(1, (K - 1) // 2 + 1):
            if math.gcd(m, K) != 1:
                continue
            f = m / K
            p = ((1 + b - MU * a) - (1 + b) * math.cos(2 * math.pi * f)) / den
            q = (1 - b) * math.sin(2 * math.pi * f) / den
            best = max(best, p * (1 - p) - q * q - q / math.tan(math.pi / K))
    return best


def gap_top(rho, K):
    tops = []
    for KK in (K, K + 1):
        g = lambda d: lag_one(tun(*unrolled_to_tuning(KAPPA, rho, d)), math.cos(2 * math.pi / KK))
        ds = [i * 0.001 for i in range(600)]
        for d0, d1 in zip(ds, ds[1:]):
            if g(d0) > 0 >= g(d1):
                tops.append(brentq(g, d0, d1, xtol=1e-12))
                break
    return min(tops)


def find(K, n_rho=80):
    k = 2 * K + 1
    lo, hi = 1 / (K + 1), 1 / K
    best = (-math.inf, None)
    for i in range(1, n_rho):
        rho = lo + (hi - lo) * i / n_rho
        dt = gap_top(rho, K)
        for j in range(80, 100):
            a, b = unrolled_to_tuning(KAPPA, rho, dt * j / 100)
            m = mixed_margin(a, b, k)
            if m > best[0]:
                best = (m, (a, b))
    obj = lambda z: -min(mixed_margin(z[0], z[1], k), -max_pure_mode(z[0], z[1]))
    res = minimize(obj, best[1], method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=400))
    a, b = res.x
    A, B = F(a).limit_denominator(1000), F(b).limit_denominator(1000)
    t = Tuning(f"fan tongue 2/{k} in gap I_{K}", A, B, KAPPA)
    mr, ar = mixed(t, k), omega_aug_report(t)
    rho, delta = tuning_to_unrolled(t)
    return dict(K_gap=K, k=k, alpha=str(A), beta=str(B), kappa=str(KAPPA), mixed_feasible=mr["feasible"],
                mixed_margin=mr["margin"], mixed_support=mr.get("support"), in_omega_aug=ar["member"], K0=ar["K0"],
                rho=rho, delta=delta)


if __name__ == "__main__":
    out = [find(K) for K in (3, 4)]
    for r in out:
        print(r)
    with open(os.path.join(os.path.dirname(__file__), "..", "results", "fan_tunings.json"), "w") as fh:
        json.dump(out, fh, indent=1)
