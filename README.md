# hb1d-cycles: evidence for the one-dimensional consequence of the Fotopoulos et al. conjecture

**Do heavy-ball cycles exist on one-dimensional smooth strongly convex functions at parameters outside the
augmented root-of-unity region Ω_aug? Can Pasechnyuk-Vilensky's multi-harmonic fan orbits be realized by a
scalar function?**

Within the tunings and cycle lengths it examines, this repository finds **no counterexample**: the answer to
both questions is no in every case checked. The results for the named tunings (§2.2) are exact and rest on two
kinds of stored evidence (§2.5). Every "no cycle" verdict is a certificate (a rational Farkas certificate, or an
interval-arithmetic certificate for the mixed condition and for Ω_aug) that stand-alone checkers outside the
package verify; every cycle found is stored as rational points that `hb1d.exact.check_cycle` and the tests
re-verify exactly. The phase diagram (§2.3) adds grid evidence whose cell classification is numerical; there,
exact 1D searches cover only the lengths stated. Nothing here proves the conjecture: one 1D cycle outside Ω_aug
would refute it, and no number of cases without one can prove it.

## 1. The question

Heavy ball, x_{k+1} = x_k − α f′(x_k) + β(x_k − x_{k−1}), on F_{μ,L}(ℝ) (μ-strongly convex, L-smooth, L = 1,
κ = L/μ).

* **Fotopoulos, Helcig, Deck, Alsbach**, *Non-Acceleration of the Heavy-Ball Method in 1D*, arXiv 2610.25104
  (18 Oct 2026). They prove Ω_aug ⊆ Ω_cyc^(1), so every tuning in Ω_aug (the union of the m-jump root-of-unity
  regions Ω_{K,m}) has a 1D cycle, and Ω_aug ⊊ Ω_Cycle (example P₁ in ℝ⁴). They **conjecture
  Ω_aug = Ω_cyc^(2)**. Since Ω_cyc^(1) ⊆ Ω_cyc^(2), the conjecture implies that **every 1D cycle lies in Ω_aug**.
* **Pasechnyuk-Vilensky**, *Two Answers on Heavy-Ball Dynamics*, arXiv 2610.09912 (7 Oct 2026). In the
  first-level gaps of the phase diagram he finds fan tongues where no single-harmonic orbit exists but
  multi-harmonic periodic orbits do (Prop. 13: κ = 100, αL = 70471/31800, β = 6169/9900, a 13-periodic orbit in
  ℝ⁵ with harmonics 2/13 and 3/13). He states as open: "the realization of the fan orbits … by scalar functions
  is the corresponding question in dimension one".

The two questions meet: as Pasechnyuk-Vilensky notes after his Thm A, his single-harmonic region
{Θ(2πj/k) ≤ π/k} is exactly Ω_aug (the one-line reason is in [docs/method.md §4](docs/method.md), checked here in
SymPy). So the fan tongues lie outside Ω_aug, and a 1D fan orbit would refute the conjecture.

## 2. Results

Tags: **[exact]** = verified by exact rational arithmetic (`fractions.Fraction`) or, where cosines of 2π/K occur,
by interval arithmetic with outward rounding (mpmath, 256 bits); **[numerical]** = floating-point LP/MILP
verdicts. Every "no cycle" below is the strongest form: it covers cycles with distinct points, with repeated
points and with slopes equal to μ or L, and all orbits whose period divides K ([docs/method.md §1](docs/method.md)).

### 2.1 Gate: the gap tunings are outside Ω_aug [exact]

`hb1d.regions.omega_aug_report` implements Fotopoulos et al.'s Lemma 3.1 and their Appendix B uniform bound. It
reproduces their Table 1 (13 rows, 9 digits), the Appendix B interval and K₀ = 10 at P₁ exactly, and Example 3.2
(in Ω_{5,2}, not in Ω_rou). It certifies that Prop. 13, conv-gap and the two fan tunings below lie **outside
Ω_aug**: every length K ≥ K₀ is excluded by the uniform bound and every primitive mode with K < K₀ has a
certified negative margin (Prop. 13: K₀ = 28, largest margin −0.00155 at (K, m) = (5, 1)).

### 2.2 Main table

| tuning (κ = L/μ) | αL | β | in Ω_aug | 1D result | how |
|---|---|---|---|---|---|
| **Prop. 13** (κ = 100), 13-cycle in ℝ⁵ | 70471/31800 | 6169/9900 | no | **no 1D cycle of any length K ≤ 25** [exact] | exact 1D search for K = 3..14; mixed condition certified infeasible for K = 15..25 (and for every K ≤ 14 except 13) |
| conv-gap (κ = 100), converges | 111/50 | 31/50 | no | no cycle of length K ≤ 40 in any dimension [exact] | mixed condition certified infeasible for K ≤ 40; the exact 1D search for K = 3..12 agrees |
| **FAN7** (κ = 100), 2/7 fan tongue, gap I₃ | 34/13 | 312/665 | no | **no 1D cycle for K ≤ 14** [exact] | exact 1D search for K = 3..14 |
| **FAN9** (κ = 100), 2/9 fan tongue, gap I₄ | 2201/986 | 199/320 | no | **no 1D cycle for K ≤ 16** [exact] | exact 1D search for K = 3..14; mixed condition certified infeasible at K = 15, 16 |
| **P₁** of Fotopoulos et al. (κ = 10⁴/531), cycle in ℝ⁴ | 32829/12500 | 4879/10000 | no | **no 1D cycle for K ≤ 13** [exact] | exact 1D search for K = 10, 13; mixed condition certified infeasible for every other K ≤ 13 |
| Polyak (κ = 100) | 400/121 | 81/121 | yes (3-tongue) | 1D cycles for every K = 3..8 [exact] | rational cycles, re-verified by `check_cycle` |
| Example 3.2 of Fotopoulos et al. (κ = 5/4) | 18/5 | 17/20 | yes (Ω_{5,2}) | 1D 5-cycle [exact]; none for K = 3, 4, 6 [exact] | rational cycle; Farkas certificates for K = 3, 4, 6 |

FAN7 and FAN9 are found by `scripts/find_fan_tunings.py`. At both, the mixed condition of Fotopoulos et al.
holds strictly at the fan length with modes 1/k and 2/k, so a cycle exists in ℝ⁴ (their Prop. 5.1). They are
fan-tongue points in Pasechnyuk-Vilensky's own sense too: his program of Thm fan-div(a), including the
minimizer constraint, is strictly feasible there with harmonics 1/k and 2/k, and the dips lie in (1/4, 1/3) and
(1/5, 1/4), i.e. in the gaps I₃ and I₄ [numerical, margins 1.3·10⁻³ and 1.9·10⁻³]. They are explicit
rational points with two-harmonic higher-dimensional cycles and **no** 1D cycle of length ≤ 14.

At Prop. 13 the mixed condition (a necessary condition for a K-cycle in **every** dimension) is certified
feasible for K ≤ 60 only at K ∈ {13, 26, 35, 39, 44, 48, 52, 53, 56, 57} and certified infeasible at every
other K ≤ 60 [exact]. So the **first open length for a 1D cycle at Prop. 13 is K = 26**.

Raw records (nodes, pruned nodes, seconds, solver versions) are in `results/exact_search.jsonl`; the table is
regenerated by `scripts/summarize.py` (`results/summary.md`). Every exclusion is stored as a checkable
certificate (§2.5).

### 2.3 Phase diagram in one dimension, κ = 100 [numerical classification, exact 1D searches]

![1D phase diagram at kappa = 100](figures/phase_diagram_1d_k100.png)

The grid has 1,740 tunings in Pasechnyuk-Vilensky's unrolled coordinates, over the gaps I₃, I₄, I₅ **and the
lower part of the tongues just above them**: 60 columns in the rotation number ρ times 29 heights
h = δ/δ_top(ρ). Of these heights, 25 lie in the gaps (h < 1, finer near the tops) and 4 above the gap tops
(h = 1.03, 1.06, 1.10, 1.15: the tongue cells drawn above the black line). The run (`scripts/phase_diagram.py`)
takes about 10 minutes on 4 cores. Each cell is a rational tuning (denominators ≤ 10⁵):

| cell class | cells | 1D result |
|---|---|---|
| in Ω_aug (the K-tongues, K = 3..6, above the gaps; membership rechecked rigorously in every cell) | 240 | a 1D cycle in **every** cell, of exactly the predicted minimal period K [exact] |
| `mixed_feasible`: outside Ω_aug, mixed condition feasible at some K ≤ 10 [numerical] | 34 | the exact 1D search found **no cycle** at the lengths that passed the filter (36 searches, K = 7, 9, 10) [exact]; other lengths K ≤ 10 are excluded only by the numerical filter |
| `mixed_feasible_long`: outside Ω_aug, mixed condition feasible only at some 10 < K ≤ 30 [numerical] (e.g. around Prop. 13, K = 13) | 39 | not searched on the grid |
| `nocycle`: mixed condition infeasible for all K ≤ 30 [numerical] | 1,427 | no cycle of length ≤ 30 in any dimension [numerical] |
| `undecided`: a mixed-condition margin within the tolerance 10⁻⁸, or an LP solver failure | 0 | — |

So on this grid a 1D cycle was found in every Ω_aug cell (as Fotopoulos et al.'s Thm 4.1 guarantees), and **no
1D cycle was found outside Ω_aug within the lengths examined**: K ∈ {7, 9, 10} in the 34 `mixed_feasible` cells,
no 1D search in the 39 cells whose first mixed-feasible length exceeds 10, and only lengths K ≤ 30 (numerically,
in any dimension) in the other 1,427 cells. The grid says nothing about longer cycles in those cells.

The mixed-condition classes are floating-point verdicts. Which lengths get an exact 1D search, and the existence
of higher-dimensional cycles in the `mixed_feasible` cells (Fotopoulos et al. Prop. 5.1), rest on that numerical
filter and are not certified cell by cell; the named tunings FAN7, FAN9 and Prop. 13 do have rigorous
mixed-condition certificates (§2.5). Ω_aug membership was rechecked rigorously in all 1,740 cells, and the 1D
searches inside the cells are exact.

For each cell outside Ω_aug and each K ≤ 30, `results/phase_diagram_k100.json` stores the LP record: the
verdict, both margins (σ for a feasible weight vector, τ for an infeasibility vector; by LP duality σ = −τ, which
holds to 4·10⁻¹⁵) and both solver statuses. The LP solutions themselves (the weight vectors s and y) are in
`results/phase_diagram_k100_lp_vectors.json.gz`, so the margins can be recomputed without the solver (a test does
this for a sample of cells), and the HiGHS version (1.12.0, bundled with SciPy 1.18.1) is recorded. All 42,000
LP records (84,000 LP solves, two per record) reported success, and the smallest margin was 9.6·10⁻⁸, ten times
the tolerance, so no verdict was undecided.

The class name `mixed_feasible` does not claim that a cell lies in a fan tongue in Pasechnyuk-Vilensky's sense.
His fan program (Thm fan-div(a)) was nevertheless evaluated in every such cell (field `pv_fan_program`) and is
strictly feasible in all 34 at one of the searched lengths [numerical; best margins 5·10⁻⁵ to 9·10⁻³]. The
named tunings (circles) are separate points with their own results in §2.2; the neighbouring cells are
different tunings. Pasechnyuk-Vilensky's `dmivilensky/Heavy-Ball` already classifies these gaps in unrestricted
dimension; the new part here is the one-dimensional column.

### 2.4 Validation against known results (all in `tests/`)

| known result | source | check |
|---|---|---|
| Table 1, Appendix B interval, sin²_min, K₀ = 10 at P₁ | Fotopoulos et al. | reproduced exactly |
| Example 3.2 in Ω_{5,2} \ Ω_rou; p₂(1−p₂)−q₂² = (17215√5−36285)/10368; min h = 63/500 | Fotopoulos et al. | SymPy-exact |
| Thm 4.1 projected 5-cycle at Example 3.2 | Fotopoulos et al. | slopes strictly inside [μ, L] (independent oracle) |
| P₁ weights s₂ = 137/500, s₃ = 363/500 strictly mixed-feasible | Fotopoulos et al. eq. (5.6) | interval arithmetic |
| Prop. 13 spectrum 149/797 at 2/13, 648/797 at 3/13 feasible; fails at conv-gap | Pasechnyuk-Vilensky | fan program with our Θ |
| GTD's h(αμ, c) = (μ/L)·𝓛(c); Θ(2πm/K) = −arg M_m + arg(1 − M_m) | Pasechnyuk-Vilensky eq. (gtd) | numerical, SymPy |
| Polyak κ = 100: 3-cycle −1 → 0 → 1 | Liu et al. 2610.10369, App. H | exact |
| Polyak κ = 100: tied 5-cycle (1, 0, −1, −δ/β, 1), δ = 1 − β | Hernanz, Thm "fails from rest" | exact (tie branch) |
| No cycle in any dimension at (q, a, b) = (1/100, 113/50, 3/5) | F. Wang, Trace231/heavyball-blog | mixed infeasible K ≤ 40 (notation assumed) |
| MILP and exact search agree on small K | `hb1d.milp` (numerical) vs `hb1d.exact` | 6 cases |

### 2.5 Stored certificates and how to check them without rerunning the searches

| file | what it proves | independent checker |
|---|---|---|
| `certificates/<tuning>_K<K>.farkas.gz` (one per exact "none" in §2.2) | no 1D K-cycle at that tuning | `scripts/verify_certificates.py` |
| `certificates/grid/*.farkas.gz` | no 1D cycle in each `mixed_feasible` cell of §2.3 at the searched K | same |
| `certificates/mixed_and_omega_aug.json` | the mixed-condition verdicts of §2.2 and Ω_aug membership of the named tunings and of all 1,740 grid cells | `scripts/verify_rigorous.py` |
| `results/phase_diagram_k100_lp_vectors.json.gz` | the LP solutions behind every numerical mixed-condition verdict of §2.3 (not a proof) | `tests/test_results.py` recomputes the margins from them |
| `results/exact_search.jsonl`, `results/phase_diagram_k100.json` (field `x`) | the cycles found (Polyak, Example 3.2, every Ω_aug cell) as exact rationals | `hb1d.exact.check_cycle`, `tests/` |

A `.farkas.gz` file lists, for every pruned node of the search tree, its prefix of the sorted order and the
support of its Farkas vector. `verify_certificates.py` uses only the Python standard library: it rebuilds each
node's rows, recomputes the Farkas vector exactly on the stored support, checks y ≥ 0 and the identity
B^T y = (∗, −1, …, −1) in rational arithmetic, and walks the tree to check that the certified prefixes cover
every sorted order. `verify_rigorous.py` uses only mpmath interval arithmetic and recomputes the Appendix B
bound exactly. Neither imports this package or an LP solver. The tests check every stored Farkas file with K ≤ 10
and that tampered files of both kinds are rejected (a dropped or altered Farkas certificate, a forged Ω_aug
witness, an edited bound or mixed-condition vector); the per-push CI also runs `verify_rigorous.py` and `verify_certificates.py` on
every file with K ≤ 12. The K = 13, 14 files (about 1 to 2 hours) are checked by a separate, manually started
workflow (`.github/workflows/full-certificates.yml`, Actions → "Full certificate check") and by
`FULL=1 ./reproduce.sh`; the last full local check is logged in `results/log_verify_certificates.txt`.

## 3. Method in one paragraph

A 1D K-cycle forces its gradients, g = Gx with a circulant rational G, and in one dimension interpolation by
F_{μ,L}(ℝ) is a set of slope conditions on consecutive sorted points. So a cycle exists iff one of (K−1)!
linear systems is feasible. We search the sorted orders by branch and bound; a branch is cut only when an
**exact rational Farkas certificate** of infeasibility has been verified, and any cycle found is rebuilt as an
exact rational vertex and checked. Long lengths are filtered first by the **mixed condition** of Fotopoulos et
al., a necessary condition for K-cycles in every dimension whose infeasibility we certify by interval arithmetic.
Details, the big-M derivation of the auxiliary MILP and the phase-diagram coordinates are in
[docs/method.md](docs/method.md).

## 4. What this does not show

* It does **not** prove the conjecture, nor its 1D consequence. The evidence is finite: named tunings, lengths
  K ≤ 14 by search (≤ 25 at Prop. 13 with the filter), and a grid with K ≤ 10.
* At Prop. 13, lengths K = 26, 35, … are open in 1D; K = 26 is out of reach of the present search
  (the tree grows by a factor of about 3 to 4 per unit of K).
* Nothing is claimed in dimension 2. Existence of a 1D cycle outside Ω_aug would refute the conjecture; its
  absence says nothing about 2D.
* The phase-diagram cell labels come from a floating-point mixed-condition classification. Ω_aug membership
  of the cells is rigorous, and the 1D searches inside the cells are exact.

## 5. Related work

* B. Goujaud, A. Taylor, A. Dieuleveut, *Provable non-accelerations of the heavy-ball method*, Math. Program.
  (2025), arXiv 2307.11291: root-of-unity cycles in ℝ², Ω_rou, the LP characterization of general cycling.
* B. Goujaud et al., *Open Problem: Two Riddles in Heavy-Ball Dynamics*, COLT 2025, arXiv 2502.19916: fractional
  sets; the 1D question.
* K. Fotopoulos, M. Helcig, K. Deck, J. Alsbach, arXiv 2610.25104, and their Lean repository
  `kostfoto/Non-Acceleration-Heavy-Ball-1D` (no computational search of Ω_cyc^(1); no commits since 17 Oct 2026).
  Their acknowledgments mention that Dieuleveut and collaborators concurrently resolved 1D non-acceleration
  (unpublished as of 22 Oct 2026).
* D. Pasechnyuk-Vilensky (arXiv's spelling; the paper itself writes Pasechniuk-Vilensky), *Two Answers on
  Heavy-Ball Dynamics*, arXiv 2610.09912, and `dmivilensky/Heavy-Ball`, which **already classifies the κ = 100
  first-level gaps without a dimension restriction**; our diagram adds the 1D question on top of it.
* K. Liu, Q. Wu, Z. Sun, M. Chen, S. Li, B. Qu, *Sharp Initialization Guarantees for the Heavy-Ball Method*,
  arXiv 2610.10369: App. F (single-harmonic 1D cycles, global non-acceleration for κ ≥ 24), App. H (attracting
  3-cycle at Polyak's tuning, κ = 100).
* M. Hernanz, *The Optimal Rate Constant of Heavy-Ball Momentum on Smooth Strongly Convex Functions* (alphaXiv,
  Oct 2026), `marcoshernanz/heavy-ball-optimal-rate`: App. D, results in one dimension, including an exact
  **two-harmonic** 1D 5-cycle with a repeated point at Polyak's tuning. That cycle is inside Ω_aug, so a
  multi-harmonic 1D cycle by itself is not new; the question here is about tunings outside Ω_aug.
* Lessard, Recht, Packard (2016); Taylor, Van Scoy, Lessard (2018); Badithela, Seiler (ACC 2019): Lyapunov and
  IQC convergence regions, the 3-cycle at Polyak's tuning.
* Control side: P. Seiler, J. Carrasco, *Construction of periodic counterexamples to the discrete-time Kalman
  conjecture*, IEEE L-CSS 5(4), 2021 (single-frequency 1D construction, the analogue of Fotopoulos et al.'s
  projection); J. Zhang, J. Carrasco, W. P. Heath, *Duality bounds for discrete-time Zames–Falb multipliers*,
  IEEE TAC 67(7), 2022, who state that the meaning of the multiple-frequency condition "is still open";
  A. Kharitenko, arXiv 2607.23599 (read at abstract level), a numerical example where every O'Shea–Zames–Falb
  test fails but the scalar Lur'e system is stable.

The literature check (channels, queries, verdicts) is logged in
[docs/literature-check-2026-10-21.md](docs/literature-check-2026-10-21.md).

## 6. Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[test]"
.venv/bin/python -m pytest -q          # 4 to 6 minutes (Apple M1)
./reproduce.sh                         # default part, about 22 minutes (M1); uses $PYTHON, else .venv/bin/python, else python3
FULL=1 ./reproduce.sh                  # also the long searches (K = 12..14), the phase diagram and the K = 13, 14 certificate check (hours)
.venv/bin/python scripts/verify_certificates.py --max-k 12 certificates/*.farkas.gz certificates/grid/*.farkas.gz   # ~9 min
.venv/bin/python scripts/verify_certificates.py certificates/*.farkas.gz certificates/grid/*.farkas.gz              # all, 1 to 2 hours (M1)
.venv/bin/python scripts/verify_rigorous.py      # interval-arithmetic certificates, ~1 min
```

Reruns replace the earlier record for the same (tuning, K) in `results/exact_search.jsonl` and rewrite its
certificate file (`scripts/run_exact.py --append` keeps the old records).

Python ≥ 3.11, NumPy, SciPy (HiGHS), mpmath and Matplotlib; the tests also need pytest and SymPy
(`pip install -e ".[test]"`). Results were computed with Python 3.13.5,
SciPy 1.18.1 (HiGHS 1.12.0) and mpmath 1.3.0 on an Apple M1.

## 7. ملخص بالعربية

**السؤال:** هل توجد دورات لطريقة Heavy Ball على دوال ذات بعد واحد عند معاملات خارج المنطقة Ω_aug؟ وهل يمكن
تحقيق مدارات «المروحة» متعددة التوافقيات عند Pasechnyuk-Vilensky بدالة ذات بعد واحد؟ وجود دورة واحدة كهذه
يدحض حدسية Fotopoulos وزملائه.

**النتيجة ضمن الحالات والأطوال المفحوصة:** لم نجد مثالاً مضاداً. النتائج الثلاث التالية للنقاط المسماة
`[exact]`، أي بحساب كسري أو بحساب فترات دقيق، وشهاداتها محفوظة في `certificates/` ويفحصها مدقق مستقل:
* معاملات Prop. 13 خارج Ω_aug، ولا توجد لها دورة في بعد واحد بأي طول K ≤ 25. الأطوال من 3 إلى 14 بالبحث
  الدقيق، وبقية الأطوال حتى 25 لأن الشرط المختلط (شرط لازم في كل الأبعاد) غير ممكن. أول طول مفتوح هو 26.
* نقطتان جديدتان داخل لساني المروحة 2/7 و2/9 عند κ = 100: لهما دورات في ℝ⁴، ولا دورة في بعد واحد بطول
  ≤ 14 و≤ 16 على الترتيب.
* النقطة P₁ عند Fotopoulos (دورة في ℝ⁴): لا دورة في بعد واحد بطول ≤ 13.

**مخطط الطور (التصنيف عددي):** مخطط في بعد واحد لفجوات المستوى الأول عند κ = 100 وللجزء الأدنى من الألسنة
فوق قممها (1,740 خلية). تصنيف الخلايا بالشرط المختلط عددي، أما عضوية Ω_aug فأُعيد فحصها بدقة في كل الخلايا،
والبحث في بعد واحد داخل الخلايا دقيق. وُجدت دورة في بعد واحد في كل خلية من خلايا Ω_aug (240 خلية). أما خارج
Ω_aug فلم توجد دورة ضمن الأطوال المفحوصة فقط: K ∈ {7, 9, 10} في الخلايا التي يكون فيها الشرط المختلط ممكناً
عند طول ≤ 10 (34 خلية)، ولا بحث في 39 خلية، واستبعاد عددي للأطوال K ≤ 30 في الباقي.

**ما لا تثبته هذه النتائج:** لا تثبت الحدسية. هي أدلة محدودة بأطوال ونقاط معينة، ولا تقول شيئاً عن بعدين.

## 8. Citation and license

To cite this repository, see `CITATION.cff`. Licensed under the Apache License 2.0 (`LICENSE`).
