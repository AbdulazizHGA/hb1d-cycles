"""Run the exact 1D cycle search for named tunings and a range of K; append one JSON line per (tuning, K)
to results/exact_search.jsonl.

    python scripts/run_exact.py prop13 3 14
    python scripts/run_exact.py ex32 5 5 --time-limit 600
"""
import argparse
import json
import os
import platform
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import scipy  # noqa: E402

from hb1d.exact import check_cycle, search, search_parallel  # noqa: E402
from hb1d.tunings import NAMED, Tuning  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tuning", help=f"one of {sorted(NAMED)} or alpha,beta,kappa as fractions, e.g. 111/50,31/50,100")
    ap.add_argument("kmin", type=int)
    ap.add_argument("kmax", type=int)
    ap.add_argument("--time-limit", type=float, default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results", "exact_search.jsonl"))
    args = ap.parse_args()

    if args.tuning in NAMED:
        t = NAMED[args.tuning]
    else:
        a, b, k = (F(s) for s in args.tuning.split(","))
        t = Tuning(args.tuning, a, b, k)

    for K in range(args.kmin, args.kmax + 1):
        if args.workers > 1:
            cycles, st = search_parallel(t, K, workers=args.workers)
        else:
            cycles, st = search(t, K, time_limit=args.time_limit)
        rec = dict(tuning=t.name, key=args.tuning, alpha=str(t.alpha), beta=str(t.beta), kappa=str(t.kappa), K=K,
                   verdict=("cycle" if cycles else ("none [exact]" if st["complete"] and st["unresolved"] == 0
                                                      else "undecided")),
                   **{k: v for k, v in st.items() if k not in ("tuning", "K")})
        if cycles:
            c = check_cycle(t, cycles[0])
            rec.update(x=[str(v) for v in cycles[0]], period=c["period"], strict_vertex=c["strict"], order=c["order"])
        rec.update(scipy=scipy.__version__, python=platform.python_version(), date=time.strftime("%Y-%m-%d %H:%M:%S"))
        print(json.dumps({k: rec[k] for k in ("key", "K", "verdict", "nodes", "pruned", "unresolved", "seconds")}),
              flush=True)
        with open(args.out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
