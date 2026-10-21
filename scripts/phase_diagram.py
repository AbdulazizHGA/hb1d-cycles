"""One-dimensional phase diagram of heavy-ball cycles in the first-level gaps at kappa = 100, together with the
lower part of the tongues just above the gap tops (heights h > 1).

Grid in Pasechnyuk-Vilensky's unrolled coordinates (rho, delta): rho = tangency rotation number on the boundary of
Omega_1, delta = distance in alpha L beyond that boundary. Each grid tuning is rounded to a rational (denominators
<= 10^5) and classified:

  omega1     Lc > 0 on [-1, 1]: converges (Pasechnyuk Thm A(i))                                    [numerical]
  aug        in Omega_aug (some pure mode K <= 40 has A - q cot(pi/K) >= 0); the exact 1D search is run at the
             smallest such K and is expected to find a cycle (Fotopoulos Thm 4.1)       [1D search: exact]
  mixed_feasible       outside Omega_aug, mixed condition feasible (or undecided) at some K <= K1D: exact 1D
                       search at each such K; Pasechnyuk-Vilensky's fan program is also evaluated at those K
                       (field pv_fan_program), but the class name does not claim it        [1D search: exact]
  mixed_feasible_long  outside Omega_aug, mixed condition feasible only at lengths in (K1D, KMIX]: 1D not searched
  nocycle              outside Omega_aug and Omega_1, mixed condition infeasible for every K <= KMIX: no cycle of
                       length <= KMIX in any dimension                                             [numerical]
  undecided            no mixed-feasible K <= KMIX, but some K has an undecided mixed verdict (margin within the
                       tolerance, or an LP solver failure)
  unresolved           some exact 1D search had unresolved nodes (cells on a region boundary, degenerate LPs)

Every cell outside Omega_aug stores, for each K <= KMIX, the mixed-condition LP record (verdict, margins sigma and
tau, solver statuses; hb1d.mixed.mixed_float_record).
  CANDIDATE  a 1D cycle outside Omega_aug (then re-checked with the rigorous omega_aug_report)

    python scripts/phase_diagram.py --nrho 48 --ndelta 36 --workers 4
"""
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction as F
from multiprocessing import get_context

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hb1d.exact import check_cycle, search  # noqa: E402
from hb1d.geometry import in_omega1, unrolled_to_tuning  # noqa: E402
from hb1d.mixed import MIXED_TOL, fan_program_margin, mixed_float_record  # noqa: E402
from hb1d.regions import omega_aug_report  # noqa: E402
from hb1d.tunings import Tuning  # noqa: E402
from scripts.run_exact import write_certs  # noqa: E402

GRID_CERTS = os.path.join(os.path.dirname(__file__), "..", "certificates", "grid")

KAPPA = F(100)
KMIX, K1D, KAUG = 30, 10, 40


def pure_mode_witness(a, b, mu=0.01):
    for K in range(3, KAUG + 1):
        for m in range(1, (K - 1) // 2 + 1):
            if math.gcd(m, K) != 1:
                continue
            den = (1 - mu) * a
            f = m / K
            p = ((1 + b - mu * a) - (1 + b) * math.cos(2 * math.pi * f)) / den
            q = (1 - b) * math.sin(2 * math.pi * f) / den
            if p * (1 - p) - q * q - q / math.tan(math.pi / K) >= 0:
                return K, m
    return None


def gap_top(rho):
    """delta at which the vertical line rho of gap I_K (K = floor(1/rho)) enters the K- or (K+1)-tongue."""
    from scipy.optimize import brentq
    from hb1d.geometry import lag_one
    K = int(1 / rho)
    tops = []
    for KK in (K, K + 1):
        g = lambda d: lag_one(Tuning("g", *map(F, unrolled_to_tuning(KAPPA, rho, d)), KAPPA), math.cos(2 * math.pi / KK))
        ds = [i * 0.001 for i in range(600)]
        for x0, x1 in zip(ds, ds[1:]):
            if g(x0) > 0 >= g(x1):
                tops.append(brentq(g, x0, x1, xtol=1e-12))
                break
    return min(tops)


# h = 1 is the tongue boundary itself (degenerate cycles, not exactifiable in floating point), so it is skipped.
H_LEVELS = [0.05 + 0.08 * i for i in range(10)] + [0.86 + 0.01 * i for i in range(14)] + [0.995, 1.03, 1.06, 1.10, 1.15]


def classify(args):
    rho, delta = args
    a, b = unrolled_to_tuning(KAPPA, rho, delta)
    t = Tuning(f"grid rho={rho:.5f} delta={delta:.5f}", F(a).limit_denominator(10 ** 5), F(b).limit_denominator(10 ** 5),
               KAPPA)
    a, b = float(t.alpha), float(t.beta)
    out = dict(rho=rho, delta=delta, alpha=str(t.alpha), beta=str(t.beta))
    if not (abs(b) < 1 and 0 < a < 2 * (1 + b)):
        out["cls"] = "outside_cv"
        return out
    if in_omega1(t):
        out["cls"] = "omega1"
        return out
    w = pure_mode_witness(a, b)
    if w is not None:
        K, m = w
        cyc, st = search(t, K)
        out.update(cls="aug" if cyc else "unresolved", witness=[K, m], found_1d=bool(cyc), nodes=st["nodes"],
                   period=check_cycle(t, cyc[0])["period"] if cyc else None,
                   x=[str(v) for v in cyc[0]] if cyc else None)
        return out
    recs = {K: mixed_float_record(t, K) for K in range(3, KMIX + 1)}
    out["mixed"] = recs
    feas = [K for K, r in recs.items() if r["v"] == 1]
    undec = [K for K, r in recs.items() if r["v"] == 0]
    out["mixed_feasible"] = feas
    out["mixed_undecided"] = undec
    if not feas:
        out["cls"] = "undecided" if undec else "nocycle"
        return out
    short = sorted(K for K in feas + undec if K <= K1D)
    if not short:
        out["cls"] = "mixed_feasible_long"
        return out
    out["pv_fan_program"] = {K: fan_program_margin(t, K)[0] for K in short}
    res = {}
    for K in short:
        certs = []
        cyc, st = search(t, K, record=certs)
        if cyc:
            rep = omega_aug_report(t)
            out.update(cls="CANDIDATE" if rep["member"] is False else "aug_boundary", K_cycle=K,
                       x=[str(v) for v in cyc[0]], omega_aug_member=rep["member"])
            return out
        res[K] = dict(complete=st["complete"], unresolved=st["unresolved"], nodes=st["nodes"])
        if st["complete"] and st["unresolved"] == 0:
            key = f"grid_rho{rho:.6f}_delta{delta:.6f}"
            res[K]["certificate"] = write_certs(GRID_CERTS, key, t, K, certs,
                                                dict(nodes=st["nodes"], pruned=st["pruned"], date="grid"))
    ok = all(v["complete"] and v["unresolved"] == 0 for v in res.values())
    out.update(cls="mixed_feasible" if ok else "unresolved", exact_1d=res)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rho", type=float, nargs=2, default=[1 / 6, 1 / 3])
    ap.add_argument("--delta", type=float, nargs=2, default=[0.002, 0.10])
    ap.add_argument("--nrho", type=int, default=48)
    ap.add_argument("--ndelta", type=int, default=36)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--normalized", action="store_true",
                    help="vertical coordinate h = delta / delta_top(rho) on the fixed levels H_LEVELS")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "results", "phase_diagram_k100.json"))
    args = ap.parse_args()
    r0, r1 = args.rho
    d0, d1 = args.delta
    rhos = [r0 + (r1 - r0) * (i + 0.5) / args.nrho for i in range(args.nrho)]
    if args.normalized:
        pts = []
        for r in rhos:
            try:
                dt = gap_top(r)
            except ValueError:          # rho at a phase-locked point 1/K: the gap closes
                continue
            pts += [(r, h * dt) for h in H_LEVELS]
    else:
        pts = [(r, d0 + (d1 - d0) * j / (args.ndelta - 1)) for r in rhos for j in range(args.ndelta)]
    t0 = time.time()
    with get_context("spawn").Pool(args.workers) as pool:
        cells = []
        for n, c in enumerate(pool.imap(classify, pts, chunksize=4), 1):
            cells.append(c)
            if c["cls"] == "CANDIDATE":
                print("CANDIDATE:", c, flush=True)
            if n % 100 == 0:
                print(f"{n}/{len(pts)} cells, {time.time() - t0:.0f}s", flush=True)
    import scipy
    meta = dict(mixed_tol=MIXED_TOL, normalized=args.normalized, H_LEVELS=H_LEVELS if args.normalized else None, kappa=str(KAPPA), KMIX=KMIX, K1D=K1D, KAUG=KAUG, nrho=args.nrho, ndelta=args.ndelta, rho=args.rho,
                delta=args.delta, seconds=time.time() - t0, scipy=scipy.__version__,
                date=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(args.out, "w") as fh:
        json.dump(dict(meta=meta, cells=cells), fh)
    from collections import Counter
    print(Counter(c["cls"] for c in cells))


if __name__ == "__main__":
    main()
