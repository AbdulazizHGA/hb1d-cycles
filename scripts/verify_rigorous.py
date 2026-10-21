"""Independent checker for certificates/mixed_and_omega_aug.json. Needs only the standard library and mpmath
(interval arithmetic with outward rounding, 256 bits); it does not import hb1d.

Mixed condition (Fotopoulos et al. eq. (5.4)), with f = m/K, den = (1 - mu) alpha (L = 1),
    p = ((1 + beta - mu alpha) - (1 + beta) cos 2 pi f)/den,  q = (1 - beta) sin 2 pi f/den,  A = p(1 - p) - q^2,
    phi(f, D) = A (cos 2 pi f D - 1) + q sin 2 pi f D.
  infeasible at K: the stored y >= 0 has sum_D y_D phi(m/K, D) > 0 for every m = 1..K-1 (theorem of alternatives);
  feasible at K:   the stored s >= 0, not all 0, has sum_m s_m phi(m/K, D) < 0 for every D = 1..K-1.

Omega_aug (Fotopoulos et al. Lemma 3.1 and Appendix B), for a tuning in Omega_cv:
  R = 1/(16 c^2 sin_min^2) with c = (1 - beta)/den and sin_min^2 = 1 - max(lo^2, hi^2), [lo, hi] = {u in [-1, 1] :
  p(u) in [0, 1]}, recomputed exactly; outside: cot(pi/K0)^2 > R and A - q cot(pi/K) < 0 for every primitive mode
  1 <= m <= (K-1)/2 with K < K0; member: the witness mode has A - q cot(pi/K) >= 0 (certified > 0).

    python scripts/verify_rigorous.py
"""
import json
import os
import sys
from fractions import Fraction as F
from math import gcd

from mpmath import iv

iv.prec = 256


def I(x):
    x = F(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def pq(alpha, beta, mu, m, K):
    den = (1 - mu) * alpha
    ang = 2 * iv.pi * m / K
    p = ((1 + beta - mu * alpha) - (1 + beta) * iv.cos(ang)) / den
    q = (1 - beta) * iv.sin(ang) / den
    return p, q, ang


def phi(alpha, beta, mu, m, K, D):
    p, q, ang = pq(alpha, beta, mu, m, K)
    A = p * (1 - p) - q * q
    return A * (iv.cos(ang * D) - 1) + q * iv.sin(ang * D)


def check_mixed(t, K, e):
    a, b, mu = I(t["alpha"]), I(t["beta"]), 1 / I(t["kappa"])
    if e["feasible"] is False:
        y = [F(v) for v in e["y"]]
        if any(v < 0 for v in y):
            return False
        return all(sum((I(y[D - 1]) * phi(a, b, mu, m, K, D) for D in range(1, K) if y[D - 1]), iv.mpf(0)).a > 0
                   for m in range(1, K))
    if e["feasible"] is True:
        s = [F(v) for v in e["s"]]
        if any(v < 0 for v in s) or not any(s):
            return False
        return all(sum((I(s[m - 1]) * phi(a, b, mu, m, K, D) for m in range(1, K) if s[m - 1]), iv.mpf(0)).b < 0
                   for D in range(1, K))
    return True          # undecided entries claim nothing


def mode_value(a, b, mu, K, m):
    p, q, _ = pq(a, b, mu, m, K)
    return p * (1 - p) - q * q - q * iv.cot(iv.pi / K)


def check_aug(e):
    alpha, beta, kappa = F(e["alpha"]), F(e["beta"]), F(e["kappa"])
    mu = 1 / kappa
    if not (abs(beta) < 1 and 0 < alpha < 2 * (1 + beta)):
        return e["member"] is False
    a, b, m_ = I(alpha), I(beta), I(mu)
    if e["member"] is True:
        K, m = e["witness"]
        return gcd(m, K) == 1 and mode_value(a, b, m_, K, m).a > 0
    if e["member"] is not False:
        return True
    den = (1 - mu) * alpha
    P0, P1 = (1 + beta - mu * alpha) / den, (1 + beta) / den
    lo, hi = max(F(-1), (P0 - 1) / P1), min(F(1), P0 / P1)
    smin = 1 - max(lo * lo, hi * hi)
    c = (1 - beta) / den
    R = 1 / (16 * c * c * smin)
    K0 = e["K0"]
    if F(e["R"]) != R or not smin > 0 or not (iv.cot(iv.pi / K0) ** 2 - I(R)).a > 0:
        return False
    return all(mode_value(a, b, m_, K, m).b < 0 for K in range(3, K0) for m in range(1, (K - 1) // 2 + 1)
               if gcd(m, K) == 1)


def main(path):
    d = json.load(open(path))
    assert d["format"] == "hb1d-rigorous-v1"
    bad = []
    for key, v in d["named"].items():
        if not check_aug(v["omega_aug"]):
            bad.append((key, "omega_aug"))
        for K, e in v["mixed"].items():
            if not check_mixed(v["omega_aug"], int(K), e):
                bad.append((key, "mixed", K))
        n_inf = sum(e["feasible"] is False for e in v["mixed"].values())
        n_fea = sum(e["feasible"] is True for e in v["mixed"].values())
        print(f"{key}: Omega_aug member = {v['omega_aug']['member']}; mixed infeasible at {n_inf} lengths, "
              f"feasible at {n_fea}", flush=True)
    gbad = [i for i, e in enumerate(d["grid"]) if not check_aug(e)]
    print(f"grid: {len(d['grid']) - len(gbad)}/{len(d['grid'])} Omega_aug verdicts verified")
    bad += [("grid", i) for i in gbad]
    print("OK" if not bad else f"FAILED: {bad[:10]}")
    return not bad


if __name__ == "__main__":
    here = os.path.dirname(__file__)
    sys.exit(0 if main(sys.argv[1] if len(sys.argv) > 1 else
                       os.path.join(here, "..", "certificates", "mixed_and_omega_aug.json")) else 1)
