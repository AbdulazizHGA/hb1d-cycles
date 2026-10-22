# Literature / novelty check: 1D heavy-ball cycles outside Omega_aug; 1D realization of fan orbits

First run 2026-10-21; last re-checked 2026-10-22 (see the end of this file). This is a snapshot: re-run the
queries below before any new release.

**How this log was made.** The searches were run on 2026-10-21 by an AI research assistant (Claude, Anthropic)
working for the repository author, using web search, direct page fetches, the Semantic Scholar and OpenAlex APIs
and the GitHub API; the author reviewed the result. Shorter re-checks were run on 2026-10-21 just before
publication and on 2026-10-22 (end of this file). Every [read] item was read in primary form; items marked
[2nd-hand] or [abs] were not read in full.

**Limits.** A search that finds nothing does not show that no prior work exists. Google Scholar could not be
queried and the arXiv API was rate-limited on 2026-10-21 (arXiv listing pages were used instead); indexes lag
behind arXiv by days to weeks; unpublished work (e.g. the concurrent work by Dieuleveut et al. mentioned in
Fotopoulos et al.'s acknowledgments) cannot be found this way.

## Direct links and queries (to re-run the check)

| what | link |
|---|---|
| Fotopoulos, Helcig, Deck, Alsbach, arXiv 2610.25104 | https://arxiv.org/abs/2610.25104 |
| Pasechnyuk-Vilensky, arXiv 2610.09912 | https://arxiv.org/abs/2610.09912 |
| Liu, Wu, Sun, Chen, Li, Qu, arXiv 2610.10369 | https://arxiv.org/abs/2610.10369 |
| Goujaud, Taylor, Dieuleveut, arXiv 2307.11291 | https://arxiv.org/abs/2307.11291 |
| Goujaud et al., COLT 2025 open problem, arXiv 2502.19916 | https://arxiv.org/abs/2502.19916 |
| Kharitenko, arXiv 2607.23599 | https://arxiv.org/abs/2607.23599 |
| Seiler, Carrasco, arXiv 2009.02468 | https://arxiv.org/abs/2009.02468 |
| Zhang, Carrasco, Heath, arXiv 2008.11975 | https://arxiv.org/abs/2008.11975 |
| Heath, Carrasco, Zhang, arXiv 2109.00744 | https://arxiv.org/abs/2109.00744 |
| Fotopoulos et al. Lean repository | https://github.com/kostfoto/Non-Acceleration-Heavy-Ball-1D |
| Pasechnyuk-Vilensky code | https://github.com/dmivilensky/Heavy-Ball |
| Hernanz repository (paper on alphaXiv) | https://github.com/marcoshernanz/heavy-ball-optimal-rate |
| F. Wang's notes | https://github.com/Trace231/heavyball-blog |
| ETH ODI blog post (7 Jul 2026) | https://odi.inf.ethz.ch/blog/heavy-ball-one-dimension/ |
| citing papers, 2610.25104 | https://api.semanticscholar.org/graph/v1/paper/arXiv:2610.25104/citations?fields=title,externalIds |
| citing papers, 2610.09912 | https://api.semanticscholar.org/graph/v1/paper/arXiv:2610.09912/citations?fields=title,externalIds |
| newest heavy-ball papers on arXiv | https://arxiv.org/search/?query=heavy-ball&searchtype=all&order=-announced_date_first |
| new commits, Fotopoulos et al. repository | `gh api repos/kostfoto/Non-Acceleration-Heavy-Ball-1D/commits` |
| new commits, Hernanz repository | `gh api repos/marcoshernanz/heavy-ball-optimal-rate/commits` |

The remaining queries (arXiv searches for "Lur'e", "Zames-Falb", "Kalman conjecture", author search
"Dieuleveut", GitHub repository and code searches) are listed with their results in the sections below.

**Scope.** The question "do HB cycles exist on 1D functions in F_{mu,L}(R) at (alpha, beta) outside Omega_aug
(Fotopoulos–Helcig–Deck–Alsbach)?", equivalently "can Pasechnyuk-Vilensky's multi-harmonic fan orbits be
realized by scalar functions?". Focus: kappa = 100, level-1 gaps.

**Legend.** **[read]**: the primary text was read (LaTeX source, PDF text or repository file). **[abs]**: the
abstract page only. **[2nd-hand]**: someone else's statement or a fetch-tool summary only.

---

## BOTTOM LINE

**Nothing I found answers the question, in either form.** No paper, preprint, repo, or blog post (as of
2026-10-21) exhibits a 1D cycle at a tuning outside Omega_aug. None realizes a fan orbit of
Pasechnyuk-Vilensky (Thm fan-div) by a scalar function, and none proves that this is impossible.

- **The two anchors still stand.** Fotopoulos et al. 2610.25104 state the conjecture
  Omega_aug = Omega_cyc^(2) (a strict subset of Omega_cyc) [read]. Pasechnyuk-Vilensky 2610.09912 state the
  1D realization of fan orbits as the open question in dimension one [read]. Neither paper has a v2, and
  no citing paper turned up.
- **Scope warning 1: 1D multi-harmonic cycles already exist in the literature, but inside the cycling
  region.** Hernanz (Thm "Polyak's tuning fails from rest", App. D) gives an exact 1D 5-cycle
  (1, 0, -1, -delta/beta, 1) for Polyak's tuning, kappa >= 100 [read]. My DFT at kappa = 100 gives
  |c_1| = 0.562 and |c_2| = 0.058, so the cycle has **two nonzero harmonics**. It also uses a repeated point.
  Polyak's tuning at kappa = 100 already has single-harmonic 3-cycles (Liu et al. App. H; our MILP).
  It is also inside Omega_aug, verified numerically with Fotopoulos' angle condition
  |arg M_1| + arg(1 - M_1) <= pi/K (quoted in the blog and in PV §7). At mu = 1, L = 100, gamma = 4/121,
  beta = 81/121 the left side is 0.458, below pi/3 for K = 3. It is also 0.399 < pi/4 for K = 4 and
  0.419 < pi/5 for K = 5. **So we cannot claim "first multi-harmonic 1D cycle". Any claim
  must be scoped to tunings outside Omega_aug, or inside the Omega_1 gaps.** A statement
  such as "Hernanz has no fan orbits" is correct in letter (his cycle is not a fan orbit and not outside
  Omega_aug), but it should say that Hernanz does have a two-harmonic 1D cycle, at Polyak's tuning.
- **Scope warning 2: the generic control analogue has been answered negatively.** Kharitenko
  2607.23599 (Jul 2026) reports a numerical disproof of Carrasco's conjecture "no OZF multiplier =>
  not absolutely stable" for discrete-time scalar slope-restricted Lur'e systems [2nd-hand: fetch-tool
  summary of the HTML]. So in general, the failure of every Zames–Falb multiplier does **not** imply that a
  destabilizing scalar nonlinearity exists. This supports the motivation that 1D realizability of fan
  obstructions is a real question. It does not answer the heavy-ball case, which has a specific plant and
  exact interpolation (Pasechnyuk-Vilensky's framing). On the realizability of multi-frequency duality
  certificates, Zhang–Carrasco–Heath (TAC 2022) say "an interpretation of the multiple frequency condition
  is still open" [read]. I found no later work that closes this.
- **Unconfirmed concurrent work.** Fotopoulos' acknowledgment says Dieuleveut and collaborators
  independently resolved **1D non-acceleration**, communicated privately [read]. That is not our question.
  On arXiv there is no 2026 Dieuleveut, Goujaud or Taylor paper on 1D heavy-ball [read: author search].
  Hernanz reports "no such preprint as of October 15, 2026" [2nd-hand]. **This is still a risk:**
  unpublished work by that group could contain 1D cycle constructions beyond projected rotations.

---

## (a) Citation channels

| Channel | Query | Result | Verdict |
|---|---|---|---|
| Semantic Scholar | `/paper/arXiv:2610.25104/citations` | `{"offset":0,"data":[]}`, a real empty answer (not a 429). The metadata call returned 429. | Nothing prior: 0 indexed citers |
| Semantic Scholar | `/paper/arXiv:2610.09912/citations` | `data: []`. Metadata: citationCount 0, referenceCount 0. S2 has not parsed its references, so it does not appear as a citer of 2502.19916. | Nothing prior |
| Semantic Scholar | `/paper/arXiv:2610.10369/citations` | `data: []` | Nothing prior |
| Semantic Scholar | `/paper/arXiv:2307.11291/citations` | 39 citers. The 2026 heavy-ball ones are 2610.25104, 2610.16034 (Lee–Zummo–Wu, SGD with momentum time series), 2610.10369, 2609.39114 (Muon), 2609.09743 (He–Zhang, randomized schedules), 2609.08656 (Ma–Zhang, convex lower bound) and a few reviews. | Nothing on 1D cycles outside Omega_aug |
| Semantic Scholar | `/paper/arXiv:2502.19916/citations` | 2610.25104, 2610.10369, 2509.18346 (Parimi et al.) | Nothing prior |
| OpenAlex | title search | 2610.09912 and 2610.10369 indexed with 0 citations. 2610.25104 not indexed. | Nothing prior |
| Google Scholar | title query for 2610.25104 | Page returned "cannot complete the operation", so no data. | Not reachable |

**Method-adjacent, title only:** another 2307.11291 citer is **2412.11330**, Ranjan, Park, Gualandi,
Lodi, Stellato, "Exact Verification of First-Order Methods via Mixed-Integer Linear Programming" (Dec 2024).
Our approach is also a MILP. Read it before citing; not characterized here.

Hernanz's paper (alphaXiv only) cites 2610.09912, 2610.10369 and the ETH blog [read]. No index lists it as a
citer.

## (b) arXiv, Oct 2026 and the Dieuleveut / Goujaud / Taylor group

| Query (arxiv.org/search, newest first) | Result | Verdict |
|---|---|---|
| `heavy-ball` (all fields) | After 2610.25104 (18 Oct) the only heavy-ball paper is 2610.27990 (Suvroikova–Dvinskikh–Dorn–Gasnikov, momentum and Polyak–Ruppert averaging, 20 Oct). Others in Oct: 2610.24421 (distributed HB), 2610.23301 (heavy tails), 2610.10369, 2610.09912. | Nothing prior |
| `Lur'e` (all) | 2026 items: 2609.14336, 2607.23599, 2608.06915 (Gramlich–Scherer–Ebenbauer, SGD with momentum via IQC), 2605.08488, and others not relevant. Nothing after 18 Oct. | See (e) |
| `Zames-Falb` | 2610.23206 (Abdelsamed–Gramlich–Ebenbauer–Scherer, transient synthesis, 16 Oct), 2607.23599, 2608.06915, 2606.14471, 2511.12826, 2510.21490, 2503.05875 | See (e) |
| `Kalman conjecture` | The only control items are 2109.00744, 2009.02468, 1704.02484 and 2207.14579. Nothing from 2026. | See (e) |
| `heavy ball cycle`, `momentum cycle strongly convex` | First returned 429; the second returned no results | Inconclusive |
| arXiv API (export.arxiv.org) | "Rate exceeded" every time | Not usable today |
| Author search `Dieuleveut`, 2025–2026 | 2610.06432, 2607.08753 (Finding Simple Proofs for First-Order Optimization, with Goujaud and Taylor), 2606.29403, 2605.31594, 2605.22929 (Prox-TMM), 2602.03329, 2601.07021, 2510.25458. **No 1D heavy-ball paper.** [read: listing] | Nothing prior; risk remains (see bottom line) |
| adrientaylor.github.io/publications | Heavy-ball entries are only 2502.19916 and 2307.11291 (Math. Program. Ser. B, 2025) [read: page summary] | Nothing prior |
| WebSearch (extended): "heavy-ball one-dimensional cycles non-acceleration arXiv October 2026" and "Dieuleveut Goujaud Taylor heavy-ball one-dimensional non-acceleration 2026" | Only 2610.25104, 2610.09912, 2610.10369, 2307.11291, 2502.19916 | Nothing prior |
| WebSearch (extended): "heavy-ball fan orbits / multi-harmonic periodic orbits scalar one-dimensional 2026" | Off-topic results only | Nothing prior |

Version check [abs pages]: 2610.25104 v1 only (18 Oct 2026, math.OC and cs.LG). 2610.09912 v1 only
(7 Oct 2026, 63 pp.; the author is spelled **Pasechnyuk-Vilensky** on arXiv). 2610.10369 v1 only
(7 Oct 2026).

Ruled out: **2609.08656** (Ma–Zhang) is a last-iterate lower bound for HB on smooth convex (not strongly
convex) functions with schedules [abs]. It does not discuss cycles or Omega regions.

### Key papers, as read from source

**Liu, Wu, Sun, Chen, Li, Qu: "Sharp Initialization Guarantees for the Heavy-Ball Method", arXiv
2610.10369** (AISTATS-2027 preprint style). Full author list: Kang Liu, Quanen Wu, Zhan Sun, Mengxiao Chen,
Suyan Li, Bohao Qu. Affiliations: XJTU, PolyU, HUST, NUS, A*STAR. [read: LaTeX source in
stage1/data/src/2610.10369]
Appendix letters: A Related Work, B Energy, C Startup, D Smooth Cycles, E Three-Cycle Constants,
**F Global Non-Acceleration**, G Necessary Cycle Size, **H Cycles from Zero Velocity**, I Numerical
Protocol.

- **App. F (Theorem "Global non-acceleration").** For kappa >= 24:
  (kappa-12)/(kappa+12) <= rho_* <= (kappa-1)/(kappa+1), also on C^infty and Hessian-Hölder subclasses.
  Every tuning with 12(1-beta) <= gamma L < 2(1+beta) has a C^infty 1D objective with a locally attracting
  K-cycle. The construction states that it builds on Fotopoulos' projected cycles. It uses
  x_j = cos(j theta - theta/4), theta = 2 pi/K (fixed phase theta/4, K >= 4 is the smallest with
  a_K <= 2h/3), and gradients from the cycle identity. All secants lie in
  a_K +/- 4d cos^2(pi/2K) cos(pi/K). The orbit is smoothed with flat transitions, then rescaled for the
  Hölder bound. **These are single-harmonic (single-jump, m = 1) cycles.** My reading: by the iff in
  Hernanz's projected-rotation lemma, they lie in the (closure of the) ordered root-of-unity region, so
  inside Omega_aug. Liu does not state this.
- **App. G (Prop. "Necessary cycle scale").** A work identity gives a lower bound on the diameter of any
  nonconstant cycle when the Hessian is Hölder. It rules out nonconstant cycles on quadratics. It is
  dimension-free.
- **App. H (Prop. 10, "Zero-velocity cycling").** Polyak's tuning with mu = 1, L = 100, gamma = 4/121,
  beta = 81/121. There is a C^infty f with 1 <= f'' <= 100 and the locally attracting 3-cycle
  -1 -> 0 -> 1. Every rest start x_{-1} = x_0 in [37/50, 19/25] converges to it (rational enclosures).
  Every 3-cycle is single-harmonic, so this is consistent with Omega_aug.
- Main text Sec. 4 (Thm "Sharp attraction scale", Cor. "Explicit radius bounds") gives 1D
  single-harmonic cycles x_j = R cos(j theta - phi), theta = 2 pi m/K, then embedded in R^n.
- **Verdict: related, differs.** All its cycles are single-harmonic 1D cycles, and none lies outside
  Omega_aug. No fan or multi-harmonic orbits appear.

**Fotopoulos, Helcig, Deck, Alsbach: "Non-Acceleration of the Heavy-Ball Method in 1D", arXiv
2610.25104** [read: source]
- Omega_aug = union over K >= 3 and m in U_K of Omega_{K,m} (m-jump roots-of-unity cycles).
  Thm 4.1: Omega_{K,m}° ⊆ Omega^(1)_{K-cycle}, via the phase
  phi_* = -(1/2) arg((1-M_m) M_m) - pi/(2K).
- Cor.: Omega_aug ⊆ Omega_cyc^(1), so HB does not accelerate on F(R).
- Sec. 5 gives a higher-dimensional example with Omega_aug ⊊ Omega_cyc, and **Conjecture
  (2D-restricted cycling): Omega_aug = Omega_cyc^(2) ⊊ Omega_cyc.**
- Sec. 6 gives the Lur'e connection to Seiler–Carrasco, single-frequency only.
- Acknowledgments: Chung-En Tsai independently found a 4D counterexample similar to Sec. 5. Andrey
  Kharitenko pointed out the Lur'e connection. Dieuleveut et al. concurrently resolved 1D non-acceleration
  (private communication).
- No computational search over Omega_cyc^(1) appears.
- **Verdict: source of the question; does not answer it.**

**Pasechnyuk-Vilensky: "Two Answers on Heavy-Ball Dynamics", arXiv 2610.09912** [read: source]
- Conclusion: an open problem on the completeness of level 1. It ends with "the realization of the fan
  orbits of Theorem fan-div by scalar functions is the corresponding question in dimension one."
- Intro: multi-harmonic obstructions correspond to the duality bounds of Zhang–Carrasco–Heath.
- In Badithela–Seiler (ACC 2019, pp. 4081–4085), kappa = 9+4 sqrt 5 is the threshold "at which a scalar
  function has a three-cycle" for Polyak's tuning [2nd-hand via PV; the web search confirmed the paper
  builds a limit-cycle objective but not the value 9+4 sqrt 5].
- Repo dmivilensky/Heavy-Ball (2 commits, 7 Oct 2026) [read: README, figures/gap_data.py,
  part2/gap_examples.py]. `gap_data.py` and `fig_level1.py` classify the **level-1 gap at kappa = 100**,
  unrolled along the boundary of Omega_1, into: K-tongue, two-lag certificate, feasible fan program with
  period pK+1, other single-harmonic tongue, fan-lag LP certificate, undecided. **These classifications
  are dimension-free; there is no 1D realizability check.**
- Exact tunings in gap_examples.py, with kappa = 100:
  - conv-gap: alpha L = 111/50, beta = 31/50.
  - Prop. 13: alpha L = 70471/31800, beta = 6169/9900.
- **Verdict: source of the question; does not answer it.** Our kappa = 100 phase diagram overlaps
  his gap_100 figure in parameter range. Ours differs in asking 1D realizability (point-ordering MILP
  plus exact LP/Farkas certificates), where his uses dimension-free fan programs and lag certificates.
  The write-up must say this explicitly.

## (c) GitHub

| Target / query | Result | Verdict |
|---|---|---|
| `gh api repos/kostfoto/Non-Acceleration-Heavy-Ball-1D/commits` | 5 commits, all on 2026-10-17 (cbfb071 initial … 2705c16 "authors" by helcig). Contents: Lean file, manifest, CI workflow, .gitignore. **Nothing new.** | Nothing prior: Lean formalization of the proved part only |
| `gh api repos/marcoshernanz/heavy-ball-optimal-rate/commits` | 22 commits, 2026-10-15 21:51 to 2026-10-16 15:42. The last is b836df7, a README BibTeX fix. **No commits since 16 Oct.** | See below |
| Hernanz paper/main.tex, App. D "Results in One Dimension" [read] | D.1: projected-rotations lemma. For single-harmonic x_t = cos(omega t + psi), a phase exists iff Phi_K >= 0, and the visit order matches GTD's conjecture (checked for K = 3..39). D.2: Thm, constant 50/3 (c_1(kappa) <= 100/(6+25/kappa)). D.3: Prop. "cycle method stops at 16"; it notes that multi-jump cycles (PV §7.4) have Phi_K >= 0 for some real K, so they do not beat the bound either. **D.4: Thm "Polyak's tuning fails from rest in 1D", kappa >= 100.** f' is piecewise linear and the period-5 orbit 1, 0, -1, -delta/beta, 1 has a repeated point. **My DFT at kappa = 100: |c_0| = 0.101, |c_1| = 0.562, |c_2| = 0.058, so it is two-harmonic.** He also says shapes with a free middle point appear feasible for kappa above about 45.76 (numerical), with exact instances at kappa = 49 and 64. Open questions: whether c_1(kappa) -> 8 in 1D, which "would require" 1D functions slowing every tuning as much as spirals; whether Polyak fails from rest in 1D for all kappa >= 9+4 sqrt 5. | **Related, differs:** a non-projected, two-harmonic 1D cycle exists, but only at Polyak's tuning, which is inside the cycling region. No fan orbits and no tunings outside Omega_aug. Its open "c_1 -> 8 in 1D" question is adjacent: multi-frequency 1D functions would be one route. |
| Repo search: "heavy ball cycle", "heavy-ball cycling", "heavy ball one-dimensional", "Kalman conjecture", "Lur'e periodic", "Pasechniuk", "Omega_aug", "root-of-unity cycle", "PEP heavy ball" | Only kostfoto and marcoshernanz are relevant. Also found dmivilensky/Heavy-Ball (PV code) and Kanyooo/HeavyBall (Liu Kang, one zip, 4 Oct). | Nothing prior |
| Repo search "heavy-ball" (updated after 2026-09-01) | **Trace231/heavyball-blog** (Feiming Wang, notes on the **second** riddle, 9–10 Sep 2026) [read: notes.md excerpts]. At normalized (q, a, b) = (1/100, 113/50, 3/5) the two-point Lyapunov family fails, yet no nonconstant cycle exists in any dimension (Fourier multimode cycle exclusion), with global convergence in a box. Others found: heavyball-research.github.io (ML lab blog, unrelated) and HomebrewML/HeavyBall (optimizer library, unrelated). | **Related, differs:** it proves convergence and no cycles at one kappa = 100 point, in all dimensions. Use it as a consistency check: our 1D diagram must show no cycles there (map (a, b) to (alpha L, beta) first; the notation is not checked). |
| Code search: "heavy ball cycle interpolation", "Omega_aug", "fan orbit heavy", "heavy-ball cycle milp", "heavy ball Farkas cycle" | "Omega_aug" only matched unrelated code. The others returned nothing. | Nothing prior |

## (d) ETH ODI blog

- https://odi.inf.ethz.ch/blog/heavy-ball-one-dimension/ (7 Jul 2026, "Cycling and Non-Acceleration of
  the Heavy-Ball Method in One Dimension", same four authors) [read: fetch-tool summary]. It covers the
  projected rotated m-jump cycles, the characterization |arg M_m| + arg(1 - M_m) <= pi/K, Omega_aug being
  strictly larger than Omega_rou (example in Omega_{5,2}), fractional regions, and a Lean 4 core. It
  mentions Dieuleveut et al. by personal communication. **No conjecture, no multi-harmonic or fan orbits,
  and no addendum.**
- Blog index: there are only two later posts. "Open Problems as Active Learning in the Age of AI"
  (15 Sep 2026) lists "Cycling behavior of Heavy-Ball" as solved and verified, with no new results.
  "Scaling Almost Anything?" (22 Sep 2026) is unrelated.
- **Verdict: nothing prior beyond 2610.25104.**

## (e) Control theory: Kalman conjecture, periodic orbits of scalar Lur'e systems, ZF multipliers

| Citation | Status | What it says (relevant to our question) | Verdict |
|---|---|---|---|
| P. Seiler, J. Carrasco, "Construction of periodic counterexamples to the discrete-time Kalman conjecture", IEEE Control Systems Letters 5(4):1291–1296, 2021, DOI 10.1109/LCSYS.2020.3033443, arXiv 2009.02468. Also an ACC 2021 version, DOI 10.23919/ACC50511.2021.9483270. | [read: PDF text] | Builds a monotone or slope-restricted, generally **non-odd** nonlinearity giving a T-periodic solution whenever a **single-frequency** phase condition holds at omega = 2 pi alpha/beta (period beta or 2 beta). Interpolation is by a multi-valued monotone function; the "more than one" phrase refers to non-uniqueness of the interpolating function. There is no multi-harmonic construction. | Related, differs: single harmonic. This is the analogue of Fotopoulos' projection. |
| W.P. Heath, J. Carrasco, M. de la Sen, "Second-order counterexamples to the discrete-time Kalman conjecture", Automatica 60:140–144, 2015. Also an ECC 2015 version, pp. 981–985. | [2nd-hand: search snippets and reference lists] | Second-order discrete-time counterexamples with stable periodic solutions | Related, differs |
| S. Wang, J. Carrasco, W.P. Heath, "Phase limitations of Zames–Falb multipliers", IEEE TAC 63(4):947–959, 2018, arXiv 1704.02484 | [2nd-hand: cited in the papers I read] | Frequency-interval phase limitations (continuous time) | Related, differs |
| J. Zhang, J. Carrasco, W.P. Heath, "Duality bounds for discrete-time Zames–Falb multipliers", IEEE TAC 67(7):3521–3528, 2022, arXiv 2008.11975 | [read: v2 PDF text] | Separation and duality gives single-frequency closed-form limits and **multi-frequency limits via LP**. Conjecture 2: no LTI ZF multiplier => not absolutely stable. Conclusion: [40] (Seiler–Carrasco) connected the **single-frequency** result to periodic behaviour, "however, an interpretation of the multiple frequency condition is still open". | **Key: states our control-side question as open (2020/2022).** |
| W.P. Heath, J. Carrasco, J. Zhang, "Phase limitations of multipliers at harmonics", IEEE TAC 69(1):566–573, 2024, DOI 10.1109/TAC.2023.3271855, arXiv 2109.00744 | [read: v1 PDF text] | **Continuous-time.** A phase condition at omega and a harmonic of it (example: omega and 2 omega). It shows a third-order plant with delay is a Kalman counterexample, with the periodic response only observed in simulation. The conclusion says no multiplier does not necessarily mean not absolutely stable, and "the conjecture remains open for both continuous-time and discrete-time systems". | Related, differs: continuous time, no construction of a multi-harmonic nonlinearity |
| L. Su, P. Seiler, J. Carrasco, S.Z. Khong, "On the necessity and sufficiency of discrete-time O'Shea–Zames–Falb multipliers", Automatica 2023 (S0005109823000225), arXiv 2112.07456 | [abs] | Splits the necessity conjecture into three steps and names the bottlenecks. PV cites it as [SSCK]. | Related, differs |
| **A. Kharitenko, "Existence of stable Lur'e systems for which the O'Shea–Zames–Falb stability test fails", arXiv 2607.23599 (26 Jul 2026, v1, math.OC)** | [abs plus a fetch-tool summary of the HTML; 2nd-hand for details] | **Discrete-time**, slope-restricted (mainly monotone), oddness not assumed. A full-block, copositivity-based multiplier certifies stability where no OZF multiplier exists. This is called a "(numerical)" disproof of Carrasco's necessity conjecture: rational certificates at N = 4 and slope bound k = 2, with bisection for N = 2..6. Their "kappa" is the slope bound, not a condition number. It cites Seiler–Carrasco and ZCH, but not heavy-ball or Fotopoulos (although Fotopoulos thanks Kharitenko). It constructs no periodic signal. | **Closest related.** For a generic scalar Lur'e plant, failure of every ZF multiplier does not imply a destabilizing scalar nonlinearity. So the ZF/fan obstruction can be conservative in 1D. It does not settle the HB plant. **Not read in full for this log.** |
| H. Gyotoku, T. Yuno, Y. Ebihara, D. Peaucelle, S. Tarbouriech, V. Magron, "Detecting destabilizing nonlinearities in absolute stability analysis of discrete-time feedback systems", arXiv 2503.05875 (Feb 2025). Related: arXiv 2411.14339 (dual LMIs for static OZF). | [abs] | The dual LMI with a rank-one condition extracts a destabilizing nonlinearity and a **nonzero equilibrium**, not a periodic orbit | Related, differs |
| K. Tong, C. Grussler, M.S. Chong, "Periodic fixed-points and their algebraic characteristics in discrete-time Lur'e feedback systems", arXiv 2609.14336 (13/15 Sep 2026) | [abs] | Circulant-matrix necessary condition (sector bound from singular values) for P-periodic points with memoryless nonlinearities, for all P. Not slope-restricted; not exact search. | Related, differs. Possibly a useful necessary-condition cross-check. |
| W.P. Heath, S. Das, J. Carrasco, "Multipliers for forced Lurye systems with slope-restricted nonlinearities", arXiv 2512.24453 (Dec 2025) | [abs] | Forced systems, power gain, period-preserving solutions | Not relevant |
| A. Badithela, P. Seiler, "Analysis of the heavy-ball algorithm using integral quadratic constraints", ACC 2019, pp. 4081–4085 | [2nd-hand] | Closed-form weighted off-by-one IQC for Polyak's tuning and a limit-cycle objective | Related, differs: inside the cycling region |

Key question: has anyone shown that a multi-harmonic periodic signal that defeats every Zames–Falb
multiplier is, or is not, realizable by a scalar slope-restricted nonlinearity?

- **No, as far as I found.**
- The positive single-frequency direction is Seiler–Carrasco 2021.
- The multi-frequency direction is explicitly open in ZCH (TAC 2022) and in HCZ (TAC 2024: "remains open").
- The closest 2026 result, Kharitenko 2607.23599, shows numerically that ZF failure does not imply
  instability for some scalar discrete-time Lur'e plant. That points towards "not always realizable"
  in general, but no specific multi-harmonic signal is analysed.
- WebSearch (extended) for "discrete-time Lurye periodic solution multiple frequencies Zames–Falb
  duality counterexample" and for "discrete-time Kalman conjecture two harmonics … 2023–2025" found no
  such work.

---

## To-do flagged before publishing

1. Read Kharitenko 2607.23599 in full. Is his plant or certificate structure related to the HB plant
   G(z) = -gamma z/((1-z)(beta-z))?
2. Re-word every novelty claim: "1D cycles **outside Omega_aug** / 1D realization of **fan** orbits", not
   "1D multi-harmonic cycles". Cite Hernanz App. D (two-harmonic 1D 5-cycle at Polyak's tuning) and Liu
   App. F/H.
3. Cite PV's gap_100 classification (dmivilensky/Heavy-Ball) and state how our 1D diagram differs.
4. Add a consistency point: Feiming Wang's (q, a, b) = (1/100, 113/50, 3/5) after checking the notation.
5. Re-run the S2 and arXiv API checks just before publishing; both were rate-limited today. Check
   arXiv daily listings for 21 Oct and later.

## Status of the to-do items (2026-10-21, after the computations)

1. Kharitenko 2607.23599: still read only at abstract level ([abs]); cited only as context, not as a result.
2. Done: README scopes every claim to tunings outside Omega_aug / fan orbits and cites Hernanz App. D and Liu et
   al. App. F/H (both appear as validation tests in tests/test_exact.py).
3. Done: README states that dmivilensky/Heavy-Ball classifies the kappa = 100 gaps without a dimension
   restriction and that our diagram adds the 1D question.
4. Done: tests/test_mixed.py::test_wang_no_cycle_point_consistent, under the stated reading of the notation.
5. Done (2026-10-21, before publishing; read directly): Semantic Scholar citations of 2610.25104 and 2610.09912
   are both empty (`{"offset": 0, "data": []}`); both abs pages show only [v1]; the arXiv search for "heavy-ball"
   sorted by announcement date shows nothing after 2610.27990 (20 Oct, momentum and Polyak-Ruppert averaging,
   unrelated). No new work on 1D cycles, Omega_aug or fan orbits.

## Re-check, 2026-10-22 (before the next release; read directly)

| channel | result |
|---|---|
| Semantic Scholar citations of 2610.25104 and 2610.09912 | both `{"offset": 0, "data": []}` |
| arXiv abs pages of 2610.25104 and 2610.09912 | both v1 only |
| arXiv search "heavy-ball", newest first | one new item since 20 Oct: 2610.29723 (Zhang et al., accelerated primal-dual dynamics for bilinear saddle-point problems, 21 Oct) — unrelated |
| latest commits | kostfoto/Non-Acceleration-Heavy-Ball-1D 2705c16 (17 Oct), marcoshernanz/heavy-ball-optimal-rate b836df7 (16 Oct), dmivilensky/Heavy-Ball 2f1f6de (7 Oct): no change |

Verdict: nothing new that answers or overlaps the question. Same limits as above (no Google Scholar, unpublished
work cannot be found).

