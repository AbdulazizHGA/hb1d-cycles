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
    d = json.load(open(os.path.join(CERTS, "mixed_and_omega_aug.json")))
    assert len(d["grid"]) == 1740
    assert all((e["member"] is True) == (e["cls"] == "aug") for e in d["grid"])
