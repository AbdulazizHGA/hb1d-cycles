"""Run the exact 1D cycle search for named tunings and a range of K; append one JSON line per (tuning, K)
to results/exact_search.jsonl.

    python scripts/run_exact.py prop13 3 14
    python scripts/run_exact.py ex32 5 5 --time-limit 600
"""
import argparse
import gzip
import json
import os
import platform
import re
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import scipy  # noqa: E402

from hb1d.exact import check_cycle, search, search_parallel  # noqa: E402
from hb1d.tunings import NAMED, Tuning  # noqa: E402


def file_key(key):
    """A file-name-safe form of a tuning key. Named keys are kept; a fraction triple such as '111/50,31/50,100'
    becomes '111d50_31d50_100' ('/' -> 'd', ',' -> '_', distinct so different inputs never share a name)."""
    out = key.replace("/", "d").replace(",", "_")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", out):
        raise ValueError(f"tuning key {key!r} cannot be used in a file name")
    return out


def write_certs(directory, key, t, K, certs, rec):
    """certificates/<file_key(key)>_K<K>.farkas.gz in the format read by scripts/verify_certificates.py."""
    os.makedirs(directory, exist_ok=True)
    name = f"{file_key(key)}_K{K}.farkas.gz"
    head = dict(format="hb1d-farkas-v1", key=key, tuning=t.name, alpha=str(t.alpha), beta=str(t.beta),
                kappa=str(t.kappa), K=K, nodes=rec["nodes"], pruned=rec["pruned"], date=rec["date"])
    # mtime=0 keeps the gzip bytes reproducible
    with open(os.path.join(directory, name), "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write((json.dumps(head) + "\n").encode())
        for p, sup in sorted(certs):
            gz.write((",".join(map(str, p)) + "|" + ",".join(map(str, sup)) + "\n").encode())
    path, root = os.path.abspath(os.path.join(directory, name)), os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    return os.path.relpath(path, root) if path.startswith(root + os.sep) else path   # repo-relative when inside it


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tuning", help=f"one of {sorted(NAMED)} or alpha,beta,kappa as fractions, e.g. 111/50,31/50,100")
    ap.add_argument("kmin", type=int)
    ap.add_argument("kmax", type=int)
    ap.add_argument("--time-limit", type=float, default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results", "exact_search.jsonl"))
    ap.add_argument("--certs", default=os.path.join(os.path.dirname(__file__), "..", "certificates"),
                    help="directory for <key>_K<K>.farkas.gz (written for every exact 'none' verdict); '' disables")
    ap.add_argument("--append", action="store_true",
                    help="keep earlier records for the same (key, K); by default they are replaced")
    args = ap.parse_args()

    if args.tuning in NAMED:
        t = NAMED[args.tuning]
    else:
        a, b, k = (F(s) for s in args.tuning.split(","))
        t = Tuning(args.tuning, a, b, k)

    for K in range(args.kmin, args.kmax + 1):
        certs = [] if args.certs else None
        if args.workers > 1:
            cycles, st = search_parallel(t, K, workers=args.workers, record=certs)
        else:
            cycles, st = search(t, K, time_limit=args.time_limit, record=certs)
        rec = dict(tuning=t.name, key=args.tuning, alpha=str(t.alpha), beta=str(t.beta), kappa=str(t.kappa), K=K,
                   verdict=("cycle" if cycles else ("none [exact]" if st["complete"] and st["unresolved"] == 0
                                                      else "undecided")),
                   **{k: v for k, v in st.items() if k not in ("tuning", "K")})
        if cycles:
            c = check_cycle(t, cycles[0])
            rec.update(x=[str(v) for v in cycles[0]], period=c["period"], strict_vertex=c["strict"], order=c["order"])
        rec.update(scipy=scipy.__version__, python=platform.python_version(), date=time.strftime("%Y-%m-%d %H:%M:%S"))
        if certs is not None and rec["verdict"] == "none [exact]":
            rec["certificate"] = write_certs(args.certs, args.tuning, t, K, certs, rec)
        print(json.dumps({k: rec[k] for k in ("key", "K", "verdict", "nodes", "pruned", "unresolved", "seconds")}),
              flush=True)
        keep = []
        if os.path.exists(args.out) and not args.append:
            with open(args.out) as fh:
                keep = [l for l in fh if (json.loads(l)["key"], json.loads(l)["K"]) != (args.tuning, K)]
            with open(args.out, "w") as fh:
                fh.writelines(keep)
        with open(args.out, "a") as fh:
            fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
