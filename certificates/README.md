# certificates/

Certificates for every exact exclusion claimed in the README. Both checkers are stand-alone: they do not import
the `hb1d` package or any LP solver.

| file | proves | check with |
|---|---|---|
| `<tuning>_K<K>.farkas.gz` | no 1D heavy-ball cycle of length K (any period dividing K) at that tuning | `python scripts/verify_certificates.py <files>` (standard library only, exact rational arithmetic) |
| `grid/*.farkas.gz` | the same, for each searched length of each `mixed_feasible` cell of the phase diagram | same |
| `mixed_and_omega_aug.json` | the mixed-condition verdicts (feasible / infeasible) of the named tunings, and Ω_aug membership of the named tunings and of all 1,740 grid cells | `python scripts/verify_rigorous.py` (mpmath interval arithmetic) |

## `.farkas.gz` format

A gzip text file. The first line is a JSON header with the tuning (`alpha`, `beta`, `kappa`, as exact fractions),
the length `K` and the node counts. Every further line is one pruned node of the branch and bound:

    p0,p1,...,pr|s0,s1,...

the prefix of the sorted order (`p0 = 0` is the index of the smallest point) and the support of its Farkas vector
(row indices). The checker rebuilds the node's rows, recomputes the Farkas vector exactly on that support,
checks the Farkas identity in rational arithmetic, and checks that the certified prefixes cover every sorted
order. The method is described in [`../docs/method.md`](../docs/method.md) §2.

## `mixed_and_omega_aug.json` format

For each named tuning: the Ω_aug verdict with the data that proves it (the Appendix B bound `R` and `K0` of
Fotopoulos et al. for a non-member, or a reduced witness mode `(K, m)` for a member), and for each K the
mixed-condition verdict with its rational vector (`y` for infeasibility, `s` for feasibility). For each grid cell:
the Ω_aug verdict and its data. See [`../docs/method.md`](../docs/method.md) §3–4.

Checking time on an Apple M1: `verify_rigorous.py` about 1 minute; `verify_certificates.py` about 9 minutes for the
files with K ≤ 12 (`--max-k 12`) and 1 to 2 hours for all of them.
