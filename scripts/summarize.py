"""Combine the stored evidence for the named tunings: for each tuning and each K <= 40, the 1D verdict and its reason.
This is the data behind the main table of the README (section 2.2).

Per (tuning, K):
  cycle            an exact rational 1D cycle was found (exact search, results/exact_search.jsonl)
  none:search      exact search complete with every node pruned by a verified Farkas certificate      [exact]
  none:mixed       the mixed condition is infeasible with a rigorous certificate, so no K-cycle exists in
                   any dimension (interval arithmetic)                                                 [exact]
  open             neither

Writes results/summary.json and results/summary.md.

    python scripts/summarize.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from hb1d.mixed import mixed  # noqa: E402
from hb1d.regions import omega_aug_report  # noqa: E402
from hb1d.tunings import NAMED  # noqa: E402

HERE = os.path.dirname(__file__)
KMAX = 40


def load_exact():
    best = {}
    with open(os.path.join(HERE, "..", "results", "exact_search.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            key = (r["key"], r["K"])
            if r["verdict"] == "cycle" or (r.get("complete") and r["unresolved"] == 0):
                best[key] = r
    return best


def main():
    exact = load_exact()
    rows = []
    for key in ["prop13", "convgap", "fan7", "fan9", "p1", "ex32", "polyak100"]:
        t = NAMED[key]
        aug = omega_aug_report(t)
        per_K = {}
        for K in range(3, KMAX + 1):
            e = exact.get((key, K))
            if e and e["verdict"] == "cycle":
                per_K[K] = dict(v="cycle", period=e.get("period"))
                continue
            if e and e["verdict"] == "none [exact]":
                per_K[K] = dict(v="none:search", nodes=e["nodes"], seconds=round(e["seconds"], 1))
                continue
            m = mixed(t, K)
            if m["feasible"] is False:
                per_K[K] = dict(v="none:mixed")
            else:
                per_K[K] = dict(v="open", mixed=m["feasible"], support=m.get("support"))
        first_cycle = next((K for K in per_K if per_K[K]["v"] == "cycle"), None)
        n_none = next((K - 1 for K in per_K if per_K[K]["v"] in ("open", "cycle")), KMAX)
        rows.append(dict(key=key, name=t.name, alpha=str(t.alpha), beta=str(t.beta), kappa=str(t.kappa),
                         in_omega_aug=aug["member"], K0=aug["K0"], witnesses=aug["witnesses"][:3],
                         first_1d_cycle=first_cycle, no_1d_cycle_up_to=None if first_cycle else n_none,
                         cycles_at=[K for K in per_K if per_K[K]["v"] == "cycle"],
                         searched=[K for K in per_K if per_K[K]["v"] == "none:search"],
                         mixed_excluded=[K for K in per_K if per_K[K]["v"] == "none:mixed"],
                         open=[K for K in per_K if per_K[K]["v"] == "open"], per_K=per_K))
    with open(os.path.join(HERE, "..", "results", "summary.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    def ranges(ks):
        """Compact list of integers: [3, 4, 5, 7] -> '3-5, 7'."""
        out, ks = [], sorted(ks)
        for k in ks:
            if out and k == out[-1][1] + 1:
                out[-1][1] = k
            else:
                out.append([k, k])
        return ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b in out) or "-"

    lines = ["| tuning | alpha L | beta | kappa | in Omega_aug | 1D cycle found at K | exact 1D search, no cycle, at K | "
             "other K excluded by the certified mixed condition | no 1D cycle for every K <= | first open K |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        upto = "-" if r["first_1d_cycle"] else str(r["no_1d_cycle_up_to"])
        first_open = "-" if r["first_1d_cycle"] else (str(r["open"][0]) if r["open"] else "> " + str(KMAX))
        lines.append(f"| {r['key']} | {r['alpha']} | {r['beta']} | {r['kappa']} | {r['in_omega_aug']} | "
                     f"{ranges(r['cycles_at'])} | {ranges(r['searched'])} | {ranges(r['mixed_excluded'])} | {upto} | "
                     f"{first_open} |")
    lines.append("")
    lines.append(f"K runs over 3..{KMAX}. Keys: prop13 and convgap are the tunings of Pasechnyuk-Vilensky's "
                 "Propositions 9.4 and 9.5; p1 and ex32 are Fotopoulos et al.'s P_1 and Example 3.2.")
    md = "\n".join(lines) + "\n"
    with open(os.path.join(HERE, "..", "results", "summary.md"), "w") as fh:
        fh.write(md)
    print(md)


if __name__ == "__main__":
    main()
