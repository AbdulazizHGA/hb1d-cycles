"""Numerical MILP search for 1D heavy-ball K-cycles (HiGHS via scipy.optimize.milp). Tag: [numerical].

Variables x_0..x_{K-1} in [0, 1] (x_0 = 0 is the minimum point, sum x = 1), one binary o_ij per pair i < j
(o_ij = 1 means x_i <= x_j) and a margin t in [0, 1]. For every pair, in the orientation chosen by o_ij,

    x_j - x_i >= t,    g_j - g_i - mu (x_j - x_i) >= t,    L (x_j - x_i) - (g_j - g_i) >= t,

with big-M switching. One solve maximizing t gives three verdicts:
    t* > 0  -> a cycle with distinct points and slopes strictly inside (mu, L) (then certified exactly by exact.py);
    t* = 0  -> only degenerate cycles (ties or boundary slopes) exist;
    infeasible (status 2) -> no K-cycle of any kind;  time limit (status 1) -> no verdict.
At t = 0 the order rows read x_j >= x_i in both orientations for a tie, as they must.

Big-M. With x in [0, 1]: |g_k| <= ((1 + beta) + beta + 1)/alpha = 2 (1 + beta)/alpha, so |g_j - g_i| <=
4 (1 + beta)/alpha, |x_j - x_i| <= 1, and every switched expression plus the margin t <= 1 is bounded by
4 (1 + beta)/alpha + L + 1 =: M (used below with a factor 1.01). The order rows need M_order = 2.
Symmetry breaking: x_0 is the minimum (cyclic shift). The reflection x -> -x is not broken here.
"""
import itertools
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from .tunings import Tuning


def milp_search(t: Tuning, K: int, time_limit=120.0, M=None):
    a, b, mu, L = float(t.alpha), float(t.beta), float(t.mu), 1.0
    if M is None:
        M = 1.01 * (4 * (1 + b) / a + L + 1)
    pairs = list(itertools.combinations(range(K), 2))
    nx, no = K, len(pairs)
    n = nx + no + 1
    it = n - 1
    G = np.zeros((K, K))
    for k in range(K):
        G[k, k] += (1 + b) / a
        G[k, (k - 1) % K] -= b / a
        G[k, (k + 1) % K] -= 1 / a
    rows, lo, hi = [], [], []

    def add(xc, p, oc, tc, lb, ub):
        r = np.zeros(n)
        r[:nx] = xc
        if p is not None:
            r[nx + p] = oc
        r[it] = tc
        rows.append(r)
        lo.append(lb)
        hi.append(ub)

    for p, (i, j) in enumerate(pairs):
        d = np.zeros(nx)
        d[j], d[i] = 1, -1
        dg = G[j] - G[i]
        # o = 1 (x_i <= x_j): rows >= t ; relaxed by M (1 - o)
        add(d, p, -2.0, -1.0, -2.0, np.inf)
        add(dg - mu * d, p, -M, -1.0, -M, np.inf)
        add(L * d - dg, p, -M, -1.0, -M, np.inf)
        # o = 0 (x_j <= x_i): reversed rows >= t ; relaxed by M o
        add(-d, p, 2.0, -1.0, 0.0, np.inf)
        add(-(dg - mu * d), p, M, -1.0, 0.0, np.inf)
        add(-(L * d - dg), p, M, -1.0, 0.0, np.inf)
    e0 = np.zeros(nx)
    e0[0] = 1
    add(e0, None, 0, 0, 0.0, 0.0)
    add(np.ones(nx), None, 0, 0, 1.0, 1.0)
    c = np.zeros(n)
    c[it] = -1.0
    integrality = np.concatenate([np.zeros(nx), np.ones(no), [0]])
    t0 = time.time()
    res = milp(c=c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=integrality,
               bounds=Bounds(np.zeros(n), np.ones(n)), options={"time_limit": time_limit})
    out = dict(K=K, status=int(res.status), seconds=time.time() - t0, M=M, time_limit=time_limit)
    if res.x is not None:
        x = res.x[:nx]
        out.update(margin=float(res.x[it]), x=x.tolist(), order=sorted(range(K), key=lambda k: (x[k], k)))
    if res.status == 2:
        out["verdict"] = "infeasible"
    elif res.status == 1:
        out["verdict"] = "time limit (no verdict)" if res.x is None else "time limit (incumbent only)"
    elif res.status == 0:
        out["verdict"] = "strict cycle (t* > 0)" if out["margin"] > 1e-9 else "only degenerate (t* = 0)"
    else:
        out["verdict"] = f"status {res.status}"
    return out
