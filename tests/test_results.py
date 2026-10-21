"""The stored evidence behind every claim in the README (without rerunning the long searches)."""
import json
import os

import pytest

from hb1d.exact import check_cycle
from hb1d.mixed import mixed
from hb1d.tunings import NAMED, PROP13

RES = os.path.join(os.path.dirname(__file__), "..", "results")


def records():
    out = {}
    with open(os.path.join(RES, "exact_search.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            out.setdefault((r["key"], r["K"]), []).append(r)
    return out


CLAIMED_NONE = ([("prop13", K) for K in range(3, 15)] + [("convgap", K) for K in range(3, 13)]
                + [("fan7", K) for K in range(3, 15)] + [("fan9", K) for K in range(3, 15)]
                + [("p1", 10), ("p1", 13), ("ex32", 3), ("ex32", 4)])


@pytest.mark.parametrize("key,K", CLAIMED_NONE)
def test_claimed_exact_nonexistence_is_recorded(key, K):
    recs = records()[(key, K)]
    assert any(r["verdict"] == "none [exact]" and r["complete"] and r["unresolved"] == 0 for r in recs)
    assert not any(r["verdict"] == "cycle" for r in recs)


@pytest.mark.parametrize("key,K", [("polyak100", K) for K in range(3, 9)] + [("ex32", 5)])
def test_recorded_cycles_reverify(key, K):
    from fractions import Fraction as F
    r = next(r for r in records()[(key, K)] if r["verdict"] == "cycle")
    c = check_cycle(NAMED[key], [F(v) for v in r["x"]])
    assert c["ok"] and c["period"] == K


def test_prop13_mixed_feasible_lengths_up_to_60():
    """Rigorous side of the filter: at Prop. 13, (MIX_K) is certified infeasible for every K <= 60 outside
    {13, 26, 35, 39, 44, 48, 52, 53, 56, 57}, and certified (strictly) feasible at those ten lengths."""
    feas = {13, 26, 35, 39, 44, 48, 52, 53, 56, 57}
    for K in range(3, 61):
        assert mixed(PROP13, K)["feasible"] is (K in feas), K


def test_phase_diagram_has_no_candidate_and_tongues_cycle():
    path = os.path.join(RES, "phase_diagram_k100.json")
    if not os.path.exists(path):
        pytest.skip("phase diagram not computed")
    cells = json.load(open(path))["cells"]
    assert not [c for c in cells if c["cls"] in ("CANDIDATE", "aug_boundary")]
    aug = [c for c in cells if c["cls"] == "aug"]
    assert aug and all(c["found_1d"] for c in aug)
    fan = [c for c in cells if c["cls"] == "fan"]
    assert all(v["complete"] and v["unresolved"] == 0 for c in fan for v in c["exact_1d"].values())


def test_phase_diagram_evidence_is_stored_and_checks():
    """Every Omega_aug cell stores its exact 1D cycle (re-verified here); every fan-cell search stores a
    certificate file that the independent verifier accepts."""
    from fractions import Fraction as F
    from hb1d.tunings import Tuning
    from scripts.verify_certificates import verify_file
    cells = json.load(open(os.path.join(RES, "phase_diagram_k100.json")))["cells"]
    for c in cells:
        if c["cls"] == "aug":
            t = Tuning("cell", F(c["alpha"]), F(c["beta"]), F(100))
            chk = check_cycle(t, [F(v) for v in c["x"]])
            assert chk["ok"] and chk["period"] == c["witness"][0]
        if c["cls"] == "fan":
            for K, v in c["exact_1d"].items():
                r = verify_file(os.path.join(RES, "..", v["certificate"]))
                assert r["ok"] and r["K"] == int(K)
