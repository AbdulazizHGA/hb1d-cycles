"""Exact decision procedure for heavy-ball K-cycles on one-dimensional functions in F_{mu,L}(R).

Reduction. Along a K-periodic HB orbit x_0, ..., x_{K-1} (indices mod K) the gradients are forced,

    g_k = ((1 + beta) x_k - beta x_{k-1} - x_{k+1}) / alpha,

and (x_k, g_k) is interpolable by some f in F_{mu,L}(R) iff, after sorting the points, every consecutive pair
x_i <= x_j satisfies mu (x_j - x_i) <= g_j - g_i <= L (x_j - x_i) (Taylor-Hendrickx-Glineur; Fotopoulos et al.
Thm 2.2). Ties x_i = x_j are allowed and force g_i = g_j. The conditions are homogeneous and translation
invariant, so we normalize x_0 = min x (cyclic shift) and sum x = 1 (excludes the constant orbit).

For a fixed sorted order sigma (sigma[0] = 0) the conditions are a linear system B x >= 0, x_0 = 0, sum x = 1.
By Farkas' lemma it is infeasible iff there is y >= 0 with (B^T y)_j = -1 for every j >= 1 (then
0 <= y^T B x = -sum_{j>=1} x_j = -1). Such a y is an exact rational certificate.

Search. Branch and bound over prefixes sigma[0..r] of the sorted order. A prefix is a valid relaxation of all its
completions with the rows: consecutive prefix pairs, and (sigma[r], j) for every j not in the prefix (x_j >= x_sigma[r]
and the slope condition, which holds for all pairs in 1D). A node is pruned only when an exact Farkas certificate
has been verified in rational arithmetic. A full order whose system is feasible yields a cycle, which is then
reconstructed exactly (rational vertex) and verified with Fractions.

So "no K-cycle" from `search` is an exact statement (tag [exact]) whenever stats['unresolved'] == 0: it covers
cycles with distinct points, repeated points and boundary slopes alike, for this K (minimal period or not).
"""
import time
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog

from .tunings import Tuning


def grad_matrix(t: Tuning, K: int):
    """Rational K x K matrix G with g = G x."""
    a, b = F(t.alpha), F(t.beta)
    G = [[F(0)] * K for _ in range(K)]
    for k in range(K):
        G[k][k] += (1 + b) / a
        G[k][(k - 1) % K] += -b / a
        G[k][(k + 1) % K] += -1 / a
    return G


def pair_rows(G, mu, L, i, j):
    """Rows r with r . x >= 0 encoding x_i <= x_j and mu <= slope(i, j) <= L."""
    K = len(G)
    d = [F(0)] * K
    d[j] += 1
    d[i] -= 1
    dg = [G[j][c] - G[i][c] for c in range(K)]
    return [d,
            [dg[c] - mu * d[c] for c in range(K)],
            [L * d[c] - dg[c] for c in range(K)]]


def node_rows(G, mu, L, prefix, remaining):
    rows = []
    for u, v in zip(prefix, prefix[1:]):
        rows += pair_rows(G, mu, L, u, v)
    last = prefix[-1]
    for j in remaining:
        rows += pair_rows(G, mu, L, last, j)
    return rows


# ------------------------------------------------------------------ exact linear algebra over Q

def solve_exact(A, b):
    """One solution of A z = b over Q (free variables set to 0), or None if inconsistent. A: list of rows."""
    m, n = len(A), len(A[0]) if A else 0
    M = [list(A[r]) + [b[r]] for r in range(m)]
    piv_cols, row = [], 0
    for col in range(n):
        p = next((r for r in range(row, m) if M[r][col] != 0), None)
        if p is None:
            continue
        M[row], M[p] = M[p], M[row]
        inv = 1 / M[row][col]
        M[row] = [v * inv for v in M[row]]
        for r in range(m):
            if r != row and M[r][col] != 0:
                f = M[r][col]
                M[r] = [vr - f * vp for vr, vp in zip(M[r], M[row])]
        piv_cols.append(col)
        row += 1
        if row == m:
            break
    if any(M[r][n] != 0 and all(v == 0 for v in M[r][:n]) for r in range(m)):
        return None
    z = [F(0)] * n
    for r, col in enumerate(piv_cols):
        z[col] = M[r][n]
    return z


def verify_farkas(rows, y, K):
    """Exact check: y >= 0 and (sum_r y_r rows_r)_j = -1 for all j >= 1."""
    if any(v < 0 for v in y):
        return False
    comb = [sum((y[r] * rows[r][j] for r in range(len(rows)) if y[r]), F(0)) for j in range(K)]
    return all(comb[j] == -1 for j in range(1, K))


def farkas_certificate(rows, K, tol=1e-10):
    """Find an exact Farkas certificate of infeasibility for {rows x >= 0, x_0 = 0, sum x = 1}, or None.

    Solves the alternative LP  y >= 0, (B^T y)_j = -1 (j >= 1)  with dual simplex (vertex solution), then
    recomputes y exactly on its support."""
    B = np.array([[float(v) for v in r] for r in rows])
    res = linprog(np.zeros(len(rows)), A_eq=B[:, 1:].T, b_eq=-np.ones(K - 1), bounds=(0, None), method="highs-ds")
    if res.status != 0:
        return None, res.status
    S = [r for r in range(len(rows)) if res.x[r] > tol]
    A = [[rows[r][j] for r in S] for j in range(1, K)]
    z = solve_exact(A, [F(-1)] * (K - 1))
    if z is None:
        return None, "inconsistent"
    y = [F(0)] * len(rows)
    for r, v in zip(S, z):
        y[r] = v
    return (y if verify_farkas(rows, y, K) else None), 0


# ------------------------------------------------------------------ exact cycles

def rational_vertex(rows, K, tol=1e-9):
    """Solve {rows x >= 0, x_0 = 0, sum x = 1} in floating point (vertex), then rebuild an exact rational point
    from its active set and verify it. Returns a list of Fractions or None."""
    B = np.array([[float(v) for v in r] for r in rows])
    Aeq = np.zeros((2, K))
    Aeq[0, 0] = 1
    Aeq[1, :] = 1
    res = linprog(np.zeros(K), A_ub=-B, b_ub=np.zeros(len(rows)), A_eq=Aeq, b_eq=[0, 1], bounds=(None, None),
                  method="highs-ds")
    if res.status != 0:
        return None
    x = res.x
    act = [r for r in range(len(rows)) if abs(B[r] @ x) < tol]
    A = [rows[r] for r in act] + [[F(1)] + [F(0)] * (K - 1), [F(1)] * K]
    z = solve_exact(A, [F(0)] * len(act) + [F(0), F(1)])
    if z is not None and all(sum(r[c] * z[c] for c in range(K)) >= 0 for r in rows):
        return z
    # fall back: round the float point and test (works when the point is interior)
    for den in (10 ** 6, 10 ** 9, 10 ** 12):
        z = [F(float(v)).limit_denominator(den) for v in x]
        z[0] = F(0)
        s = sum(z)
        z = [v / s for v in z]
        if all(sum(r[c] * z[c] for c in range(K)) >= 0 for r in rows):
            return z
    return None


def minimal_period(x):
    K = len(x)
    return next(p for p in range(1, K + 1) if K % p == 0 and all(x[k] == x[(k + p) % K] for k in range(K)))


def check_cycle(t: Tuning, x):
    """Exact check that the rational points x form an HB cycle on some f in F_{mu,L}(R).
    Returns dict(ok, period, strict, gradients, order)."""
    K = len(x)
    x = [F(v) for v in x]
    G = grad_matrix(t, K)
    g = [sum(G[k][c] * x[c] for c in range(K)) for k in range(K)]
    order = sorted(range(K), key=lambda k: x[k])
    mu, L = t.mu, F(1)
    ok, strict = True, True
    for u, v in zip(order, order[1:]):
        d, dg = x[v] - x[u], g[v] - g[u]
        if d == 0:
            ok &= dg == 0
            strict = False
        else:
            ok &= mu * d <= dg <= L * d
            strict &= mu * d < dg < L * d
    nonconst = len(set(x)) > 1
    return dict(ok=ok and nonconst, period=minimal_period(x), strict=ok and strict, gradients=g, order=order)


# ------------------------------------------------------------------ branch and bound

def search(t: Tuning, K: int, max_cycles=1, time_limit=None, log_every=0, start=None, record=None):
    """Exact search for K-cycles of HB on F_{mu,L}(R).

    Returns (cycles, stats). cycles: list of exact rational point lists x (verified by check_cycle) found at full
    orders; the search stops after max_cycles. stats: nodes, pruned (each with a verified Farkas certificate),
    unresolved (nodes that could be neither pruned exactly nor shown feasible exactly; must be 0 for an exact
    verdict), complete (bool: the whole tree was explored), seconds.

    If `record` is a list, every pruned node is appended to it as (prefix, support): the prefix of the sorted order
    and the indices of the nonzero entries of its verified Farkas vector y. That is enough to rebuild y exactly
    (scripts/verify_certificates.py solves B_S^T y_S = (-1, ..., -1) on the support S)."""
    G = grad_matrix(t, K)
    mu, L = t.mu, F(1)
    stats = dict(K=K, tuning=t.name, nodes=0, pruned=0, unresolved=0, leaves_feasible=0, complete=False)
    cycles = []
    t0 = time.time()
    stack = [list(p) for p in start] if start is not None else [[0]]
    while stack:
        if time_limit is not None and time.time() - t0 > time_limit:
            stats["seconds"] = time.time() - t0
            return cycles, stats
        prefix = stack.pop()
        remaining = [j for j in range(K) if j not in prefix]
        rows = node_rows(G, mu, L, prefix, remaining)
        stats["nodes"] += 1
        if log_every and stats["nodes"] % log_every == 0:
            print(f"  K={K} nodes={stats['nodes']} pruned={stats['pruned']} depth={len(prefix)} "
                  f"t={time.time() - t0:.0f}s", flush=True)
        y, status = farkas_certificate(rows, K)
        if y is not None:
            stats["pruned"] += 1
            if record is not None:
                record.append((tuple(prefix), tuple(r for r, v in enumerate(y) if v)))
            continue
        if len(remaining) <= 1:
            full = prefix + remaining
            rows = node_rows(G, mu, L, full, [])
            x = rational_vertex(rows, K)
            if x is not None and check_cycle(t, x)["ok"]:
                stats["leaves_feasible"] += 1
                cycles.append(x)
                if len(cycles) >= max_cycles:
                    stats["seconds"] = time.time() - t0
                    return cycles, stats
            else:
                stats["unresolved"] += 1
            continue
        if status == 0:
            # the alternative LP was feasible in floating point but exactification failed: branch anyway
            stats.setdefault("inexact_branch", 0)
            stats["inexact_branch"] += 1
        for j in reversed(remaining):
            stack.append(prefix + [j])
    stats["complete"] = True
    stats["seconds"] = time.time() - t0
    return cycles, stats


def _subtree(args):
    t, K, prefix, max_cycles, rec = args
    certs = [] if rec else None
    cycles, st = search(t, K, max_cycles=max_cycles, start=[prefix], record=certs)
    return cycles, st, certs


def search_parallel(t: Tuning, K: int, workers=4, depth=3, max_cycles=1, record=None):
    """Same verdict as `search`, with the subtrees below the prefixes of length `depth` explored by a process pool.
    Prefixes of length < depth are handled serially (and may be pruned there)."""
    from multiprocessing import get_context
    G = grad_matrix(t, K)
    mu, L = t.mu, F(1)
    t0 = time.time()
    agg = dict(K=K, tuning=t.name, nodes=0, pruned=0, unresolved=0, leaves_feasible=0)
    frontier, level = [[0]], 1
    while level < depth and frontier:
        nxt = []
        for prefix in frontier:
            remaining = [j for j in range(K) if j not in prefix]
            agg["nodes"] += 1
            y, _ = farkas_certificate(node_rows(G, mu, L, prefix, remaining), K)
            if y is not None:
                agg["pruned"] += 1
                if record is not None:
                    record.append((tuple(prefix), tuple(r for r, v in enumerate(y) if v)))
            else:
                nxt += [prefix + [j] for j in remaining]
        frontier, level = nxt, level + 1
    cycles = []
    with get_context("spawn").Pool(workers) as pool:
        for cyc, st, certs in pool.imap_unordered(_subtree, [(t, K, p, max_cycles, record is not None)
                                                             for p in frontier]):
            if record is not None:
                record += certs
            for k in ("nodes", "pruned", "unresolved", "leaves_feasible"):
                agg[k] += st[k]
            cycles += cyc
            if len(cycles) >= max_cycles:
                pool.terminate()
                break
    agg["complete"] = len(cycles) < max_cycles
    agg["seconds"] = time.time() - t0
    return cycles[:max_cycles], agg
