"""Omega_aug membership (Fotopoulos et al. Lemma 3.1, Appendix B) and the identity with Pasechnyuk's Theta."""
import math
import random
from fractions import Fraction as F
from math import gcd

import pytest
import sympy as sp
from mpmath import mp

from hb1d.geometry import lag_one
from hb1d.regions import (fhda_angle, in_cv, large_K_bound, omega_aug_report, pure_mode, theta_closed, theta_def)
from hb1d.tunings import CONV_GAP, FAN7, FAN9, FHDA_EX32, FHDA_P1, POLYAK100, PROP13, Tuning

# Fotopoulos et al., Table 1 (Appendix B), point P_1: (K, m, A, q cot(pi/K), A - q cot(pi/K)), 9 printed digits.
TABLE_B = [
    (3, 1, 0.101657820, 0.102961064, -0.001303244), (4, 1, 0.205813010, 0.205922128, -0.000109119),
    (5, 1, 0.191293486, 0.269555566, -0.078262080), (5, 2, -0.041604599, 0.166594501, -0.208199100),
    (6, 1, 0.152186257, 0.308883193, -0.156696936), (7, 1, 0.114644817, 0.334312475, -0.219667658),
    (7, 2, 0.178944050, 0.416880838, -0.237936788), (7, 3, -0.095868914, 0.185529427, -0.281398340),
    (8, 1, 0.083760067, 0.351531062, -0.267770994), (8, 3, 0.012302067, 0.351531062, -0.339228995),
    (9, 1, 0.059339548, 0.363667631, -0.304328082), (9, 2, 0.205071765, 0.557171135, -0.352099370),
    (9, 4, -0.120318113, 0.193503504, -0.313821617),
]


def mid(x):
    return float(mp.mpf(x.mid))


@pytest.mark.parametrize("row", TABLE_B)
def test_table_b_reproduced(row):
    K, m, A, qc, v = row
    a, q, d = pure_mode(FHDA_P1, K, m)
    assert abs(mid(a) - A) < 1e-9 and abs(mid(q) - qc) < 1e-9 and abs(mid(d) - v) < 1e-9


def test_p1_appendix_b_bound():
    # Appendix B: u in [-56921/74395, 168555301/185987500], sin_min^2 = 6180460661049399/34591350156250000, K0 = 10
    R, lo, hi, smin = large_K_bound(FHDA_P1)
    assert (lo, hi) == (F(-56921, 74395), F(168555301, 185987500))
    assert smin == F(6180460661049399, 34591350156250000)
    r = omega_aug_report(FHDA_P1)
    assert r["K0"] == 10 and r["member"] is False


def test_example_3_2_in_omega_5_2_not_rou():
    r = omega_aug_report(FHDA_EX32)
    assert r["member"] is True and r["witnesses"] == [(5, 2)]
    assert all(s == -1 for (K, m, v, s) in r["modes"] if m == 1)          # not in any Omega_{K,1}
    # exact value p_2(1 - p_2) - q_2^2 = (17215 sqrt5 - 36285)/10368 and the strict tan(pi/5) inequality
    mu, b, g = sp.Rational(4, 5), sp.Rational(17, 20), sp.Rational(18, 5)
    c, s = sp.cos(4 * sp.pi / 5), sp.sin(4 * sp.pi / 5)
    p = ((1 + b - mu * g) - (1 + b) * c) / ((1 - mu) * g)
    q = (1 - b) * s / ((1 - mu) * g)
    A = sp.nsimplify(sp.radsimp(sp.expand(p * (1 - p) - q ** 2)))
    assert sp.simplify(A - (17215 * sp.sqrt(5) - 36285) / 10368) == 0
    assert (sp.tan(sp.pi / 5) * A - q).evalf(50) > 0
    # min over c in [-1/2, 1] of GTD's h(72/25, c) equals 63/500 > 0
    cc = sp.symbols("c")
    v = mu * g
    h = v ** 2 - 2 * (b - cc + mu * (1 - b * cc)) * v + 2 * mu * (1 - cc) * (1 + b ** 2 - 2 * b * cc)
    crit = [r for r in sp.solve(sp.diff(h, cc), cc) if -sp.Rational(1, 2) <= r <= 1]
    assert min(h.subs(cc, x) for x in crit + [-sp.Rational(1, 2), 1]) == sp.Rational(63, 500)


@pytest.mark.parametrize("t", [PROP13, CONV_GAP, FAN7, FAN9])
def test_gap_tunings_outside_omega_aug(t):
    r = omega_aug_report(t)
    assert r["in_cv"] and r["member"] is False and not r["undecided"]
    assert all(s == -1 for (*_, s) in r["modes"])


def test_polyak_in_three_tongue():
    r = omega_aug_report(POLYAK100)
    assert r["member"] is True and (3, 1) in r["witnesses"]


def test_theta_identity_and_membership_agree():
    """Pasechnyuk's Theta(2 pi m/K) equals Fotopoulos' -arg M_m + arg(1 - M_m), and Theta <= pi/K <=> Omega_{K,m}."""
    rng = random.Random(0)
    n = 0
    for _ in range(60):
        kap = F(rng.choice([2, 10, 100, 1000]))
        be = F(rng.randint(1, 999), 1000)
        t = Tuning("r", F(rng.randint(1, 1999), 1000) * (1 + be), be, kap)
        if not in_cv(t):
            continue
        for K in range(3, 16):
            for m in range(1, (K - 1) // 2 + 1):
                if gcd(m, K) != 1:
                    continue
                w = 2 * math.pi * m / K
                assert abs(theta_def(t, w) - theta_closed(t, w)) < 1e-10
                assert abs(theta_def(t, w) - fhda_angle(t, K, m)) < 1e-10
                v = mid(pure_mode(t, K, m, 64)[2])
                if abs(theta_def(t, w) - math.pi / K) > 1e-9:
                    assert (theta_def(t, w) <= math.pi / K) == (v >= 0)
                n += 1
    assert n > 1000


def test_gtd_identity():
    """Pasechnyuk-Vilensky eq. (7.1): GTD's h(alpha mu, c) = (mu/L) Lc(c)."""
    rng = random.Random(1)
    for _ in range(200):
        t = Tuning("r", F(rng.randint(1, 3000), 1000), F(rng.randint(0, 999), 1000), F(rng.choice([2, 10, 100])))
        c = rng.uniform(-1, 1)
        mu, be, v = float(t.mu), float(t.beta), float(t.alpha * t.mu)
        h = v * v - 2 * (be - c + mu * (1 - be * c)) * v + 2 * mu * (1 - c) * (1 + be * be - 2 * be * c)
        assert abs(h - mu * lag_one(t, c)) < 1e-12


def test_theta_identity_symbolic():
    """Why Theta(2 pi m/K) = -arg M_m + arg(1 - M_m): with z = e^{i omega},
    (L - mu) gamma M = -a(z)/z and (L - mu) gamma (1 - M) = b(z)/z (L = 1, x = gamma mu, y = gamma), so the angle
    from M to 1 - M equals the angle from -a(z) to b(z), which is Pasechnyuk's Theta. Exact (SymPy)."""
    z, g, b, mu = sp.symbols("z gamma beta mu")
    M = ((1 + b - mu * g) - z - b / z) / ((1 - mu) * g)
    a_poly = z ** 2 - (1 + b - g * mu) * z + b
    b_poly = z ** 2 - (1 + b - g) * z + b
    assert sp.simplify((1 - mu) * g * M + a_poly / z) == 0
    assert sp.simplify((1 - mu) * g * (1 - M) - b_poly / z) == 0
