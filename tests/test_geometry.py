"""Exact Omega_1 test and solver-failure handling of the mixed-condition routine."""
import random
from fractions import Fraction as F
from types import SimpleNamespace

import pytest

import hb1d.mixed as mixed_mod
from hb1d.geometry import in_omega1, in_omega1_sampled, lag_one, lag_one_coeffs
from hb1d.regions import in_cv
from hb1d.tunings import CONV_GAP, POLYAK100, PROP13, Tuning


def test_coefficients_match_float_formula():
    rng = random.Random(5)
    for _ in range(200):
        t = Tuning("r", F(rng.randint(1, 3000), 1000), F(rng.randint(0, 999), 1000), F(rng.choice([2, 10, 100])))
        a2, a1, a0 = lag_one_coeffs(t)
        c = rng.uniform(-1, 1)
        assert abs(float(a2) * c * c + float(a1) * c + float(a0) - lag_one(t, c)) < 1e-9


def test_exact_agrees_with_sampling_away_from_the_boundary():
    """On random tunings whose minimum of Lc over [-1, 1] is not tiny, the exact test and the old 4001-point
    sampling agree."""
    rng = random.Random(6)
    n = 0
    for _ in range(400):
        be = F(rng.randint(0, 999), 1000)
        t = Tuning("r", F(rng.randint(1, 1999), 1000) * (1 + be), be, F(rng.choice([2, 10, 100, 1000])))
        if not in_cv(t):
            continue
        a2, a1, a0 = lag_one_coeffs(t)
        cands = [F(-1), F(1)] + ([-a1 / (2 * a2)] if a2 > 0 and abs(a1) < 2 * a2 else [])
        if abs(float(min((a2 * c + a1) * c + a0 for c in cands))) < 1e-6:
            continue
        assert in_omega1(t) == in_omega1_sampled(t)
        n += 1
    assert n > 200


def test_exact_catches_a_dip_between_samples():
    """A tuning with a negative dip of Lc narrower than the sample spacing 2/4000: the exact test reports 'not in
    Omega_1', the sampled one is fooled. Built by moving beta off a boundary point of Omega_1."""
    from hb1d.geometry import boundary_point
    a, b = boundary_point(100, 0.2 + 1 / 8000)   # tangency between two sample points
    t = Tuning("near boundary", F(a).limit_denominator(10 ** 12) * (1 + F(1, 10 ** 9)),
               F(b).limit_denominator(10 ** 12), F(100))
    a2, a1, a0 = lag_one_coeffs(t)
    v = -a1 / (2 * a2)
    if (a2 * v + a1) * v + a0 >= 0:
        pytest.skip("rounding put the tuning on the convergent side; nothing to compare")
    assert not in_omega1(t)
    assert in_omega1_sampled(t)


@pytest.mark.parametrize("t,expect", [(PROP13, False), (CONV_GAP, False), (POLYAK100, False),
                                      (Tuning("small step", F(1, 10), F(1, 10), F(100)), True)])
def test_known_tunings(t, expect):
    assert in_omega1(t) is expect


def test_mixed_reports_solver_failure(monkeypatch):
    """If HiGHS does not return success, mixed() returns feasible=None with the status instead of reading .fun."""
    fail = SimpleNamespace(status=4, message="numerical difficulties (simulated)", fun=None, x=None)
    monkeypatch.setattr(mixed_mod, "linprog", lambda *a, **k: fail)
    r = mixed_mod.mixed(PROP13, 13)
    assert r["feasible"] is None and r["status"] == 4 and "failed" in r["note"]
    rec = mixed_mod.mixed_float_record(PROP13, 13)
    assert rec["v"] == 0 and rec["status_feasible_lp"] == 4
