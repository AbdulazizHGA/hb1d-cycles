"""The mixed condition (necessary for K-cycles in every dimension) and Pasechnyuk's fan program."""
from fractions import Fraction as F

import pytest

from hb1d.mixed import check_weights, fan_program_slacks, mixed
from hb1d.tunings import CONV_GAP, FAN7, FAN9, FHDA_P1, POLYAK100, PROP13


def test_p1_published_weights_strictly_feasible():
    """Fotopoulos et al. eq. (5.6): s_{2/10} = 137/500, s_{3/10} = 363/500 strictly satisfy all nine inequalities."""
    ok, worst = check_weights(FHDA_P1, 10, {2: F(137, 500), 3: F(363, 500)})
    assert ok and worst < 0


def test_prop13_pasechnyuk_weights():
    """Pasechnyuk-Vilensky Prop. 13: nu = 149/797 at 2/13 and 648/797 at 3/13 satisfies the fan program; the same
    weights fail at the convergent tuning of Prop. conv-gap."""
    nu = {2: F(149, 797), 3: F(648, 797)}
    assert min(fan_program_slacks(PROP13, 13, nu)) > 1e-4
    assert min(fan_program_slacks(CONV_GAP, 13, nu)) < 0


def test_prop13_mixed_only_at_13_below_26():
    """Rigorous: at Prop. 13 the mixed condition holds at K = 13 (support {2, 3}) and fails for every other
    K <= 25, so Prop. 13 has no cycle of length K in {3..25} \\ {13} in ANY dimension."""
    for K in range(3, 26):
        r = mixed(PROP13, K)
        if K == 13:
            assert r["feasible"] is True and set(r["support"]) <= {2, 3, 10, 11}
        else:
            assert r["feasible"] is False, K


def test_convgap_mixed_infeasible():
    for K in range(3, 31):
        assert mixed(CONV_GAP, K)["feasible"] is False


@pytest.mark.parametrize("t,k", [(FAN7, 7), (FAN9, 9)])
def test_fan_tunings_mixed_feasible(t, k):
    r = mixed(t, k)
    assert r["feasible"] is True and set(r["support"]) <= {1, 2, k - 1, k - 2}
    for K in range(3, k):
        assert mixed(t, K)["feasible"] is False


def test_polyak_mixed_feasible_at_3():
    assert mixed(POLYAK100, 3)["feasible"] is True


def test_wang_no_cycle_point_consistent():
    """Feiming Wang's notes (github.com/Trace231/heavyball-blog) report no nonconstant cycle in any dimension at
    (q, a, b) = (1/100, 113/50, 3/5); read as (mu/L, alpha L, beta) (our assumption about the notation), the
    mixed condition fails for every K <= 40 there, so there is no cycle of length <= 40 in any dimension."""
    from hb1d.tunings import Tuning
    t = Tuning("Wang", F(113, 50), F(3, 5), F(100))
    assert all(mixed(t, K)["feasible"] is False for K in range(3, 41))
