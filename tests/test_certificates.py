"""The stored certificates and the two independent verifiers (scripts/verify_certificates.py, verify_rigorous.py)."""
import glob
import gzip
import json
import os

import pytest

from hb1d.exact import search
from hb1d.tunings import CONV_GAP, PROP13
from scripts.run_exact import write_certs
from scripts.verify_certificates import verify_file
from scripts import verify_rigorous

ROOT = os.path.join(os.path.dirname(__file__), "..")
CERTS = os.path.join(ROOT, "certificates")


def make(tmp_path, t, key, K):
    rec = []
    cycles, st = search(t, K, record=rec)
    assert not cycles and st["complete"]
    write_certs(str(tmp_path), key, t, K, rec, dict(nodes=st["nodes"], pruned=st["pruned"], date="test"))
    return os.path.join(tmp_path, f"{key}_K{K}.farkas.gz")


def rewrite(path, f):
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    with gzip.open(path, "wt") as fh:
        fh.writelines(f(lines))


def test_verifier_accepts_fresh_certificates(tmp_path):
    for K in (5, 7, 8):
        assert verify_file(make(tmp_path, PROP13, "prop13", K))["ok"]


def test_verifier_rejects_missing_branch(tmp_path):
    p = make(tmp_path, CONV_GAP, "convgap", 7)
    rewrite(p, lambda ls: ls[:1] + ls[2:])           # drop one certified prefix
    r = verify_file(p)
    assert not r["ok"] and r["uncovered"] > 0


def test_verifier_rejects_wrong_support(tmp_path):
    p = make(tmp_path, CONV_GAP, "convgap", 7)

    def corrupt(ls):
        prefix, sup = ls[1].strip().split("|")
        return ls[:1] + [prefix + "|" + sup.split(",")[0] + "\n"] + ls[2:]   # keep one support row only
    rewrite(p, corrupt)
    r = verify_file(p)
    assert not r["ok"] and r["invalid"] == 1


def test_verifier_rejects_other_tuning(tmp_path):
    """Certificates for Prop. 13 do not certify Polyak's tuning, which has cycles."""
    p = make(tmp_path, PROP13, "prop13", 5)

    def swap(ls):
        h = json.loads(ls[0])
        h.update(alpha="400/121", beta="81/121")
        return [json.dumps(h) + "\n"] + ls[1:]
    rewrite(p, swap)
    assert not verify_file(p)["ok"]


def stored_records():
    with open(os.path.join(ROOT, "results", "exact_search.jsonl")) as fh:
        return [json.loads(l) for l in fh]


def test_every_exact_none_has_a_certificate_file():
    for r in stored_records():
        if r["verdict"] == "none [exact]":
            path = os.path.join(ROOT, r["certificate"])
            with gzip.open(path, "rt") as fh:
                head = json.loads(fh.readline())
                n = sum(1 for _ in fh)
            assert (head["K"], head["alpha"], head["beta"]) == (r["K"], r["alpha"], r["beta"])
            assert n == r["pruned"] == head["pruned"]


@pytest.mark.parametrize("path", sorted(p for p in glob.glob(os.path.join(CERTS, "*.farkas.gz"))
                                        if int(p.rsplit("_K", 1)[1].split(".")[0]) <= 10))
def test_stored_certificates_verify(path):
    """Full independent check of every stored certificate file with K <= 10 (the larger ones are checked by
    reproduce.sh / scripts/verify_certificates.py; results/log_verify_certificates.txt records that run)."""
    assert verify_file(path)["ok"]


def test_rigorous_certificates_named_tunings():
    d = json.load(open(os.path.join(CERTS, "mixed_and_omega_aug.json")))
    for key, v in d["named"].items():
        assert verify_rigorous.check_aug(v["omega_aug"]), key
        for K, e in v["mixed"].items():
            if int(K) <= 30 or key == "prop13":
                assert verify_rigorous.check_mixed(v["omega_aug"], int(K), e), (key, K)
    assert [int(K) for K, e in d["named"]["prop13"]["mixed"].items() if e["feasible"]] == \
        [13, 26, 35, 39, 44, 48, 52, 53, 56, 57]


def test_rigorous_grid_membership_matches_classification():
    """The rigorous Omega_aug verdict of every grid cell is decided (never None) and agrees with the class stored
    in results/phase_diagram_k100.json for the same cell (matched by rho, delta and the rational tuning)."""
    d = json.load(open(os.path.join(CERTS, "mixed_and_omega_aug.json")))
    cells = json.load(open(os.path.join(ROOT, "results", "phase_diagram_k100.json")))["cells"]
    assert len(d["grid"]) == len(cells) == 1740
    for e, c in zip(d["grid"], cells):
        assert (e["rho"], e["delta"], e["alpha"], e["beta"]) == (c["rho"], c["delta"], c["alpha"], c["beta"])
        assert e["member"] in (True, False)
        assert e["member"] == (c["cls"] == "aug")


def test_fraction_input_writes_a_checkable_certificate(tmp_path):
    """run_exact.py documents fraction input such as 111/50,31/50,100; its certificate must get a safe file name
    (no '/' subfolders) and pass the independent verifier."""
    import subprocess
    import sys
    out, certs = tmp_path / "x.jsonl", tmp_path / "certs"
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "run_exact.py"), "111/50,31/50,100", "5", "5",
                    "--out", str(out), "--certs", str(certs)], check=True, capture_output=True)
    rec = json.loads(out.read_text().splitlines()[-1])
    assert rec["verdict"] == "none [exact]"
    path = certs / "111d50_31d50_100_K5.farkas.gz"
    assert path.exists() and rec["certificate"] == str(path)
    assert verify_file(str(path))["ok"]


def test_file_key_is_safe_and_injective():
    from scripts.run_exact import file_key
    assert file_key("prop13") == "prop13"
    assert file_key("111/50,31/50,100") == "111d50_31d50_100"
    assert file_key("1/1,1") != file_key("1,1/1")
    assert os.sep not in file_key("../x")             # '/' never reaches the file name
    with pytest.raises(ValueError):
        file_key("1 /2,3")                             # anything else unsafe is refused


def test_certificate_headers_match_the_named_tunings():
    """Each stored certificate's header names the exact tuning it certifies; it must be the named tuning."""
    from hb1d.tunings import NAMED
    for r in stored_records():
        if r["verdict"] != "none [exact]":
            continue
        with gzip.open(os.path.join(ROOT, r["certificate"]), "rt") as fh:
            head = json.loads(fh.readline())
        t = NAMED[r["key"]]
        assert (head["alpha"], head["beta"], head["kappa"]) == (str(t.alpha), str(t.beta), str(t.kappa))
        assert (r["alpha"], r["beta"], r["kappa"]) == (str(t.alpha), str(t.beta), str(t.kappa))


def _rigorous():
    return json.load(open(os.path.join(CERTS, "mixed_and_omega_aug.json")))


@pytest.mark.parametrize("key", ["prop13", "convgap", "fan7", "fan9", "p1"])
@pytest.mark.parametrize("witness", [[4, 3], [5, 4], [4, -1], [5, 1], [6, 1]])
def test_rigorous_verifier_rejects_forged_membership(key, witness):
    """A tuning outside Omega_aug cannot be passed off as a member with any witness, reduced or not."""
    e = dict(_rigorous()["named"][key]["omega_aug"], member=True, witness=witness)
    assert not verify_rigorous.check_aug(e)


@pytest.mark.parametrize("field,value", [("R", "1"), ("K0", 5), ("member", False)])
def test_rigorous_verifier_rejects_edited_exclusions(field, value):
    d = _rigorous()
    if field == "member":
        e = dict(d["named"]["polyak100"]["omega_aug"], member=False, K0=40, R="156")   # a member claimed outside
    else:
        e = dict(d["named"]["prop13"]["omega_aug"], **{field: value})
    assert not verify_rigorous.check_aug(e)


def test_rigorous_verifier_rejects_edited_mixed_vectors():
    from fractions import Fraction as F
    d = _rigorous()
    t = d["named"]["prop13"]["omega_aug"]
    inf = dict(d["named"]["prop13"]["mixed"]["12"])                      # infeasible at K = 12
    n = len(inf["y"])
    assert verify_rigorous.check_mixed(t, 12, inf)
    last_lag_only = ["0"] * (n - 1) + ["1"]                              # not a certificate at this tuning
    assert not verify_rigorous.check_mixed(t, 12, dict(inf, y=last_lag_only))
    assert not verify_rigorous.check_mixed(t, 12, dict(inf, y=["-1"] + inf["y"][1:]))   # negative entry
    fea = dict(d["named"]["prop13"]["mixed"]["13"])                      # feasible at K = 13
    assert verify_rigorous.check_mixed(t, 13, fea)
    assert not verify_rigorous.check_mixed(t, 13, dict(fea, s=["0"] * len(fea["s"])))
    assert not verify_rigorous.check_mixed(t, 13, dict(fea, feasible=False, y=fea["s"]))
