# results/

The authoritative records are the first five files; the certificates that prove the exact exclusions are in
[`../certificates/`](../certificates/README.md).

| file | content | written by |
|---|---|---|
| `exact_search.jsonl` | one record per (tuning, K) of the exact 1D search: verdict (`cycle` / `none [exact]`), nodes, pruned nodes, `complete`, `unresolved`, seconds, the cycle as exact rationals when one is found, the certificate file, the Python and SciPy versions (runs made after 22 Oct 2026 also record the HiGHS version) | `scripts/run_exact.py` |
| `summary.md`, `summary.json` | the main table of the README (§2.2), built from `exact_search.jsonl` and the rigorous mixed-condition verdicts | `scripts/summarize.py` |
| `phase_diagram_k100.json` | the 1,740 grid cells of the phase diagram (README §2.3): tuning, class, the exact 1D cycle of every Ω_aug cell, the exact-search records and certificate files of every `mixed_feasible` cell, and for every cell outside Ω_aug the mixed-condition LP record for each K ≤ 30 (verdict, margins, solver statuses) | `scripts/phase_diagram.py` |
| `phase_diagram_k100_lp_vectors.json.gz` | the LP solutions (weight vectors) behind the numerical mixed-condition verdicts of the grid; numerical evidence, not a proof | `scripts/phase_diagram.py` |
| `fan_tunings.json` | how the fan-tongue tunings FAN7 and FAN9 were found (margins, unrolled coordinates) | `scripts/find_fan_tunings.py` |

**Names.** The data keys `prop13` and `convgap` (in file names and records) stand for the tunings of
Pasechnyuk-Vilensky's Propositions 9.4 and 9.5 (arXiv 2610.09912; numbers as printed in the PDF). The `tuning`
labels stored inside the records and certificate headers written on 21 Oct 2026 say "Pasechnyuk Prop. 13" and
"Pasechnyuk Prop. conv-gap", after the LaTeX labels of these propositions; they are left unchanged so that the
stored files stay byte-identical to the ones that were verified.

## logs/

Console output of the runs, kept for the record. The numbers in them are the same as in the files above.

| file | run |
|---|---|
| `log_rerun_A.txt`, `log_rerun_B.txt` | the runs (21 Oct 2026) that produced the current `exact_search.jsonl` records and the stored certificates: Prop. 9.4 (`prop13`), Prop. 9.5 (`convgap`), P₁, Example 3.2 and Polyak (A); FAN9 and FAN7 (B) |
| `log_prop13_exact.txt`, `log_convgap_exact.txt`, `log_fan7_exact.txt`, `log_fan7_exact_12.txt`, `log_fan9_exact.txt`, `log_fan9_exact_13.txt`, `log_p1_exact.txt`, `log_ex32_exact.txt`, `log_polyak100_exact.txt` | the first runs (21 Oct 2026), before certificates were stored; the reruns above reproduced every verdict and node count |
| `log_phase_diagram.txt` | the last run of `scripts/phase_diagram.py` (progress lines and the class counts) |
| `log_rigorous.txt` | `scripts/export_rigorous.py` followed by `scripts/verify_rigorous.py` |
| `log_verify_certificates.txt` | the full check of all 87 Farkas certificate files (K ≤ 14) by `scripts/verify_certificates.py` |
