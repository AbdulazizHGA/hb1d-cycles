"""Export the interval-arithmetic certificates to certificates/mixed_and_omega_aug.json:

  * mixed condition (MIX_K) for the named tunings: for each K the verdict and the rational vector that proves it
    (y with sum_D y_D phi(m/K, D) > 0 for all m when infeasible; weights s with sum_m s_m phi(m/K, D) < 0 for all D
    when strictly feasible);
  * Omega_aug membership for the named tunings and for every phase-diagram cell: K0 and the exact rationals of the
    Appendix B bound, plus the witness mode for members.

Checked independently by scripts/verify_rigorous.py.

    python scripts/export_rigorous.py
"""
import json
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from hb1d.mixed import mixed  # noqa: E402
from hb1d.regions import omega_aug_report  # noqa: E402
from hb1d.tunings import NAMED, Tuning  # noqa: E402

HERE = os.path.dirname(__file__)
KMAX = {"prop13": 60}


def aug_entry(t):
    r = omega_aug_report(t)
    e = dict(alpha=str(t.alpha), beta=str(t.beta), kappa=str(t.kappa), member=r["member"])
    if r["in_cv"]:
        e.update(K0=r["K0"], R=str(r["R"]), sin_min_sq=str(r["sin_min_sq"]))
    if r["witnesses"]:
        e["witness"] = list(r["witnesses"][0])
    if r["undecided"]:
        e["undecided"] = r["undecided"]
    return e


def main():
    out = dict(format="hb1d-rigorous-v1", named={}, grid=[])
    for key, t in NAMED.items():
        mx = {}
        for K in range(3, KMAX.get(key, 40) + 1):
            m = mixed(t, K)
            entry = dict(feasible=m["feasible"])
            if m["feasible"] is False:
                entry["y"] = [str(v) for v in m["cert"]]
            elif m["feasible"] is True:
                entry["s"] = [str(v) for v in m["cert"]]
            mx[K] = entry
        out["named"][key] = dict(omega_aug=aug_entry(t), mixed=mx)
        print(key, out["named"][key]["omega_aug"]["member"], flush=True)
    path = os.path.join(HERE, "..", "results", "phase_diagram_k100.json")
    cells = json.load(open(path))["cells"]
    for c in cells:
        t = Tuning("cell", F(c["alpha"]), F(c["beta"]), F(100))
        e = aug_entry(t)
        e.update(rho=c["rho"], delta=c["delta"], cls=c["cls"])
        out["grid"].append(e)
    agree = sum((e["member"] is True) == (e["cls"] == "aug") for e in out["grid"])
    print(f"grid: rigorous Omega_aug membership agrees with the floating classification in {agree}/{len(cells)} cells")
    with open(os.path.join(HERE, "..", "certificates", "mixed_and_omega_aug.json"), "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
