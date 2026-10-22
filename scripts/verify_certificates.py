"""Independent checker for the stored "no 1D K-cycle" certificates in certificates/*.farkas.gz.

Uses only the Python standard library (fractions, gzip, json); it does not import hb1d, LP solvers or floats.

File format (gzip text): a JSON header line {"format": "hb1d-farkas-v1", "alpha", "beta", "kappa", "K", ...},
then one line per pruned node "p0,p1,...,pr|s0,s1,..." with the prefix of the sorted order and the support S of
its Farkas vector.

For a prefix (p0 = 0, ..., pr) and the remaining indices j (increasing), the node's rows r . x >= 0 are, in this
order, for each consecutive prefix pair (u, v) and then for each pair (pr, j):

    x_v - x_u,   (g_v - g_u) - mu (x_v - x_u),   L (x_v - x_u) - (g_v - g_u),

with g = G x, G_kk = (1+beta)/alpha, G_{k,k-1} = -beta/alpha, G_{k,k+1} = -1/alpha (indices mod K), L = 1,
mu = 1/kappa. Every 1D K-cycle normalized by x_0 = min x, sum x = 1 satisfies the rows of the node of its own
sorted order, and of every prefix of it (all-pairs slope conditions hold in 1D).

Checks, all in exact rational arithmetic:
  1. certificate: solve sum_{r in S} y_r rows_r[j] = -1 for j = 1..K-1 by Gaussian elimination (free variables 0)
     and require y >= 0. Then 0 <= y . (B x) = -sum_{j>=1} x_j = -1 for any feasible x: contradiction.
  2. coverage: walking the tree of prefixes from (0,), every branch reaches a certified prefix before it reaches a
     full order. Together: no sorted order admits a cycle, so there is no 1D K-cycle (of any period dividing K).

    python scripts/verify_certificates.py certificates/*.farkas.gz
    python scripts/verify_certificates.py --max-k 12 certificates/*.farkas.gz certificates/grid/*.farkas.gz
"""
import gzip
import json
import sys
from fractions import Fraction as F


def grad_matrix(alpha, beta, K):
    G = [[F(0)] * K for _ in range(K)]
    for k in range(K):
        G[k][k] += (1 + beta) / alpha
        G[k][(k - 1) % K] -= beta / alpha
        G[k][(k + 1) % K] -= 1 / alpha
    return G


def rows_for(G, mu, prefix, K):
    rem = [j for j in range(K) if j not in prefix]
    pairs = list(zip(prefix, prefix[1:])) + [(prefix[-1], j) for j in rem]
    rows = []
    for u, v in pairs:
        d = [F(0)] * K
        d[v] += 1
        d[u] -= 1
        dg = [G[v][c] - G[u][c] for c in range(K)]
        rows += [d, [dg[c] - mu * d[c] for c in range(K)], [d[c] - dg[c] for c in range(K)]]   # L = 1
    return rows


def solve(A, b):
    """Gaussian elimination over Q; one solution with free variables 0, or None."""
    m, n = len(A), len(A[0])
    M = [A[i][:] + [b[i]] for i in range(m)]
    piv, row = [], 0
    for col in range(n):
        p = next((i for i in range(row, m) if M[i][col] != 0), None)
        if p is None:
            continue
        M[row], M[p] = M[p], M[row]
        inv = 1 / M[row][col]
        M[row] = [v * inv for v in M[row]]
        for i in range(m):
            if i != row and M[i][col] != 0:
                f = M[i][col]
                M[i] = [a - f * c for a, c in zip(M[i], M[row])]
        piv.append(col)
        row += 1
        if row == m:
            break
    for i in range(row, m):
        if M[i][n] != 0:
            return None
    z = [F(0)] * n
    for i, col in enumerate(piv):
        z[col] = M[i][n]
    return z


def check_node(G, mu, K, prefix, support):
    rows = rows_for(G, mu, prefix, K)
    if any(r >= len(rows) for r in support):
        return False
    A = [[rows[r][j] for r in support] for j in range(1, K)]
    y = solve(A, [F(-1)] * (K - 1))
    if y is None or any(v < 0 for v in y):
        return False
    comb = [sum((yv * rows[r][j] for yv, r in zip(y, support)), F(0)) for j in range(K)]
    return all(comb[j] == -1 for j in range(1, K))


def verify_file(path):
    with gzip.open(path, "rt") as fh:
        head = json.loads(fh.readline())
        assert head["format"] == "hb1d-farkas-v1"
        certs = {}
        for line in fh:
            p, s = line.strip().split("|")
            certs[tuple(int(v) for v in p.split(","))] = tuple(int(v) for v in s.split(","))
    K = head["K"]
    alpha, beta, mu = F(head["alpha"]), F(head["beta"]), 1 / F(head["kappa"])
    G = grad_matrix(alpha, beta, K)
    bad = [p for p, s in certs.items() if not check_node(G, mu, K, p, s)]
    used, uncovered, stack = set(), [], [(0,)]
    while stack:
        p = stack.pop()
        if p in certs:
            used.add(p)
            continue
        rem = [j for j in range(K) if j not in p]
        if len(rem) <= 1:
            uncovered.append(p)
            continue
        stack += [p + (j,) for j in rem]
    ok = not bad and not uncovered
    # alpha, beta, kappa are part of the result: a certificate proves "no 1D K-cycle" for exactly this tuning
    return dict(file=path, key=head.get("key"), alpha=head["alpha"], beta=head["beta"], kappa=head["kappa"], K=K,
                certificates=len(certs), invalid=len(bad), uncovered=len(uncovered), unused=len(certs) - len(used),
                ok=ok)


if __name__ == "__main__":
    args = sys.argv[1:]
    max_k = None
    if args[:1] == ["--max-k"]:                # check only files whose header K is <= max_k (used by CI)
        max_k, args = int(args[1]), args[2:]
    paths = [p for p in args if not p.split("/")[-1].startswith("._")]
    if max_k is not None:
        def header_k(p):
            with gzip.open(p, "rt") as fh:
                return json.loads(fh.readline())["K"]
        skipped = [p for p in paths if header_k(p) > max_k]
        paths = [p for p in paths if p not in skipped]
        print(f"checking {len(paths)} files with K <= {max_k}; skipping {len(skipped)} larger ones", flush=True)
    results = []
    for p in paths:                     # print each file as soon as it is checked (large files take minutes)
        results.append(verify_file(p))
        print(json.dumps(results[-1]), flush=True)
    sys.exit(0 if results and all(r["ok"] for r in results) else 1)
