"""Exact 1D cycle search: known cycles are found and certified, known non-cycling tunings are refuted with
verified Farkas certificates, and the numerical MILP agrees."""
import cmath
import math
import random
from fractions import Fraction as F

import numpy as np
import pytest

from hb1d.exact import (check_cycle, farkas_certificate, grad_matrix, node_rows, search, search_parallel,
                        verify_farkas)
from hb1d.milp import milp_search
from hb1d.mixed import phi_matrix_float
from hb1d.tunings import CONV_GAP, FAN7, FAN9, FHDA_EX32, FHDA_P1, POLYAK100, PROP13, Tuning


@pytest.mark.parametrize("K", [3, 4, 5, 6])
def test_polyak_cycles_found_and_certified(K):
    cycles, st = search(POLYAK100, K)
    assert cycles, st
    c = check_cycle(POLYAK100, cycles[0])
    assert c["ok"] and c["period"] == K


def test_polyak_three_cycle_of_lessard_et_al():
    """Lessard-Recht-Packard / Liu et al. App. H: Polyak at kappa = 100 cycles on three points. Our cycle is a
    certified rational 3-cycle on F_{mu,L}(R)."""
    (x,), _ = search(POLYAK100, 3)
    assert check_cycle(POLYAK100, x)["ok"] and len(set(x)) == 3


def test_example_3_2_five_cycle():
    """Fotopoulos et al. Example 3.2 lies in Omega_{5,2}, so by their Thm 4.1 a 1D 5-cycle exists."""
    cycles, _ = search(FHDA_EX32, 5)
    assert cycles and check_cycle(FHDA_EX32, cycles[0])["ok"]


def test_theorem_4_1_projected_cycle_oracle():
    """Independent oracle: the explicit projected cycle of Fotopoulos et al. Thm 4.1 at Example 3.2 (K = 5, m = 2)
    has all adjacent sorted slopes strictly inside [mu, L] (floating point, margin >> rounding)."""
    t, K, m = FHDA_EX32, 5, 2
    al, be, mu = float(t.alpha), float(t.beta), float(t.mu)
    w = cmath.exp(2j * math.pi * m / K)
    M = ((1 + be - mu * al) - w - be / w) / ((1 - mu) * al)
    phi = -0.5 * cmath.phase((1 - M) * M) - math.pi / (2 * K)
    y = [2 * math.cos(m * k * 2 * math.pi / K + phi) for k in range(K)]
    g = [((1 + be) * y[k] - be * y[k - 1] - y[(k + 1) % K]) / al for k in range(K)]
    o = sorted(range(K), key=lambda k: y[k])
    slopes = [(g[v] - g[u]) / (y[v] - y[u]) for u, v in zip(o, o[1:])]
    assert all(mu + 1e-6 < s < 1 - 1e-6 for s in slopes)


def test_omega_rou_sample_cycles():
    """Polyak at kappa = 100 is in Omega_{4,1} and Omega_{5,1} (Omega_rou): Fotopoulos Thm 4.1 gives 1D 4- and
    5-cycles; found above. Another Omega_rou point: kappa = 100, alpha = 3, beta = 0.6 (in the 3-tongue)."""
    from hb1d.regions import omega_aug_report
    t = Tuning("rou sample", F(3), F(3, 5), F(100))
    r = omega_aug_report(t)
    assert r["member"] and (3, 1) in r["witnesses"]
    cycles, _ = search(t, 3)
    assert cycles and check_cycle(t, cycles[0])["ok"]


@pytest.mark.parametrize("t,K", [(CONV_GAP, k) for k in range(3, 9)] + [(PROP13, k) for k in range(3, 10)]
                         + [(FAN7, 7), (FAN9, 9), (FHDA_EX32, 3), (FHDA_EX32, 4)])
def test_no_cycle_exact(t, K):
    cycles, st = search(t, K)
    assert not cycles and st["complete"] and st["unresolved"] == 0


def test_parallel_matches_serial():
    _, s1 = search(PROP13, 8)
    _, s2 = search_parallel(PROP13, 8, workers=2, depth=2)
    assert s1["nodes"] == s2["nodes"] and s2["complete"] and s2["unresolved"] == 0


def test_farkas_certificates_verify_independently():
    """Re-derive certificates on random prefixes and re-check them from scratch."""
    rng = random.Random(3)
    t, K = PROP13, 9
    G = grad_matrix(t, K)
    n = 0
    for _ in range(40):
        rest = list(range(1, K))
        rng.shuffle(rest)
        prefix = [0] + rest[:rng.randint(0, K - 2)]
        remaining = [j for j in range(K) if j not in prefix]
        rows = node_rows(G, t.mu, F(1), prefix, remaining)
        y, _ = farkas_certificate(rows, K)
        if y is not None:
            assert verify_farkas(rows, y, K)
            # the combination proves 0 <= -sum_{j>=1} x_j on every feasible point: perturbing y breaks it
            y2 = list(y)
            k = next(i for i, v in enumerate(y2) if v > 0)
            y2[k] *= 2
            assert not verify_farkas(rows, y2, K) or len(set(map(tuple, rows))) < len(rows)
            n += 1
    assert n > 10


def test_check_cycle_rejects_bad_points():
    (x,), _ = search(POLYAK100, 3)
    assert not check_cycle(POLYAK100, [x[0], x[2], x[1]])["ok"] or len(set(x)) < 3
    assert not check_cycle(POLYAK100, [F(1)] * 3)["ok"]          # constant orbit is not a cycle


@pytest.mark.parametrize("t,K,expect", [(POLYAK100, 3, True), (POLYAK100, 5, True), (FHDA_EX32, 5, True),
                                        (FHDA_EX32, 4, False), (PROP13, 7, False), (CONV_GAP, 8, False)])
def test_milp_agrees_with_exact(t, K, expect):
    r = milp_search(t, K, time_limit=120)
    assert (r["verdict"] != "infeasible") == expect
    if expect:
        assert r["margin"] > 0


@pytest.mark.parametrize("K", [3, 4, 5, 6])
def test_found_cycles_satisfy_mixed_condition(K):
    """Consistency of the two theories: the DFT weights s_m = |c_m|^2 of a 1D cycle satisfy Fotopoulos' mixed
    condition sum_m s_m phi(m/K, Delta) <= 0 (it is necessary in every dimension)."""
    (x,), _ = search(POLYAK100, K)
    xs = np.array([float(v) for v in x])
    c = np.fft.fft(xs) / K
    s = np.abs(c[1:]) ** 2
    val = phi_matrix_float(POLYAK100, K) @ s
    assert np.all(val <= 1e-12)


def test_p1_gradient_matrix_rational():
    G = grad_matrix(FHDA_P1, 5)
    assert all(sum(r) == 0 for r in G)                 # translation invariance


def test_liu_et_al_appendix_h_three_cycle():
    """Liu et al. (arXiv 2610.10369) App. H: Polyak at kappa = 100 has the 3-cycle -1 -> 0 -> 1 on a 1D function."""
    c = check_cycle(POLYAK100, [F(-1), F(0), F(1)])
    assert c["ok"] and c["period"] == 3


def test_hernanz_degenerate_five_cycle():
    """Hernanz, Thm 'Polyak's tuning fails from rest in one dimension' (kappa >= 100): the period-5 orbit
    1, 0, -1, -delta/beta, 1 with delta = 1 - beta has a repeated point; our exact checker accepts it through the
    tie branch (equal points, equal gradients) and reports it as non-strict."""
    b = POLYAK100.beta
    d = 1 - b
    c = check_cycle(POLYAK100, [F(1), F(0), F(-1), -d / b, F(1)])
    assert c["ok"] and c["period"] == 5 and not c["strict"]
    assert c["gradients"][0] == c["gradients"][4]
