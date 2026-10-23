# Method

## 0. Conventions

| | this repo | Pasechnyuk-Vilensky (2610.09912) | Fotopoulos et al. (2610.25104) |
|---|---|---|---|
| method | x_{k+1} = x_k - alpha f'(x_k) + beta (x_k - x_{k-1}) | same, step alpha | same, step gamma |
| smoothness | L = 1 | general L | general L (L = 1 in Sec. 5) |
| condition number | kappa = L/mu (Prop. 9.4: kappa = 100) | kappa = L/mu | **kappa_F = mu/L** (Prop. 9.4: kappa_F = 1/100) |
| step | alpha (= alpha L, since L = 1) | alpha L | L gamma |

`hb1d.tunings` stores every tuning as exact rationals (alpha L, beta, kappa = L/mu). Fotopoulos et al.'s
Example 3.2 (kappa_F = 4/5) is stored as kappa = 5/4, and their point P_1 (kappa_F = 531/10^4) as kappa = 10^4/531.

## 1. A 1D K-cycle is a union of linear programs

Along a K-periodic orbit x_0, ..., x_{K-1} (indices mod K) the recursion forces the gradients,

    g_k = ((1 + beta) x_k - beta x_{k-1} - x_{k+1}) / alpha,

so g = G x for a circulant rational matrix G with zero row sums. In one dimension, (x_k, g_k) is the data of
some f in F_{mu,L}(R) iff, after sorting the points, every consecutive pair x_i <= x_j satisfies

    mu (x_j - x_i) <= g_j - g_i <= L (x_j - x_i)

(Taylor-Hendrickx-Glineur interpolation; Fotopoulos et al. Thm 2.2 for the 1D form). Ties x_i = x_j are allowed
and force g_i = g_j, which is what a function needs. The conditions are homogeneous and translation invariant,
so we normalize x_0 = min x (cyclic shift) and sum x = 1 (excludes the constant orbit).

For a fixed sorted order sigma with sigma[0] = 0 this is a linear system

    B_sigma x >= 0,  x_0 = 0,  sum x = 1,

and a 1D K-cycle exists iff one of the (K-1)! systems is feasible. The rows are **closed** inequalities, so one
system covers cycles with distinct points, cycles with repeated points (Hernanz's 5-cycle is an example) and
cycles with slopes equal to mu or L. A distinction between "strict" cycles (distinct points, slopes strictly
inside (mu, L)) and "non-strict" ones is therefore unnecessary: every "no K-cycle" in this repository is the
non-strict statement, which is the strongest one.

A search at length K finds every orbit whose period divides K, so "no K-cycle" also excludes, for instance,
3-cycles when K = 6.

## 2. Exact branch and bound with Farkas certificates (`hb1d.exact`)

Nodes are prefixes sigma[0..r] of the sorted order. A prefix is relaxed by the rows of its consecutive pairs and
of the pairs (sigma[r], j) for every j not yet placed (x_j >= x_{sigma[r]} and the slope condition, which holds
for all pairs in 1D because consecutive slopes in [mu, L] average to slopes in [mu, L]).

By Farkas' lemma the node system is infeasible iff there is y >= 0 with (B^T y)_j = -1 for every j >= 1: then
0 <= y^T B x = -sum_{j >= 1} x_j = -1 on any feasible point. We find y with the dual simplex (a vertex, small
support), recompute it **exactly** on its support by Gaussian elimination over Q, and check y >= 0 and
B^T y = (*, -1, ..., -1) with `fractions.Fraction`. A node is pruned **only** after this check passes. A full
order whose system cannot be refuted is solved exactly (rational vertex from the active set) and the cycle is
re-verified by `check_cycle` with Fractions. Nodes that are neither pruned nor shown feasible exactly are counted
as `unresolved`; every verdict reported below has `unresolved = 0` and `complete = true`.

Every pruned node is stored as (prefix, support of y) in `certificates/<tuning>_K<K>.farkas.gz`.
`scripts/verify_certificates.py` recomputes y from its support by exact elimination and rechecks every
certificate, and the coverage of the tree, with the standard library only.

Tag: **[exact]** (rational arithmetic). The tree size grows by a factor of about 3 to 4 per unit of K: Prop. 9.4
needs 63,362 nodes and 469 s at K = 13 on one core of an M1, and 222,835 nodes at K = 14.

## 3. The mixed condition: a necessary condition in every dimension (`hb1d.mixed`)

Write a K-cycle in any dimension through its DFT, z_t = sum_m c_m omega_K^{mt}. Fotopoulos et al. (Sec. 3 and 5)
transform the interpolation inequalities to the convex case and sum the inequality for the pair (j, j + Delta)
over j: the function values cancel and so do all cross terms between modes, leaving

    sum_m s_m phi(m/K, Delta) <= 0   (Delta = 1, ..., K-1),   s_m proportional to |c_m|^2,          (MIX_K)

(their eq. (21); eq. (8) is the two-dimensional case), with
phi(f, Delta) = A(f)(cos(2 pi f Delta) - 1) + q(f) sin(2 pi f Delta) and A = p(1 - p) - q^2 (their eq. (11)).
So if (MIX_K) has no solution s >= 0 (not all zero), there is no K-cycle in any dimension, in particular none in
1D. We only use this direction. (Fotopoulos et al. Prop. 5.1 gives the converse in dimension 2r, r = |supp s|.
Pasechnyuk-Vilensky's fan program in Theorem 9.1(a) has an extra constraint from the minimizer; it is not used as
a filter here, but to validate the angle function and, as additional information, in the grid cells.)

By the theorem of alternatives (MIX_K) is infeasible iff some y >= 0 has sum_Delta y_Delta phi(m/K, Delta) > 0 for
every m. We take y from an LP, round it to a rational vector and check every sum with **interval arithmetic**
(mpmath `iv`, 256 bits, outward rounding); the entries cos(2 pi m Delta/K) are algebraic, not rational.
Tag: **[exact]** in the sense of a rigorous computer-assisted check, the same standard as Pasechnyuk-Vilensky's
Prop. 9.4 ("proved with exact rational data ... interval arithmetic").

This filter is what makes long lengths reachable: at Prop. 9.4, (MIX_K) is feasible for K <= 60 only at
K = 13, 26, 35, 39, 44, 48, 52, 53, 56, 57, so the exact 1D search is needed only at those K.

## 4. Omega_aug membership (`hb1d.regions`)

Fotopoulos et al. Lemma 3.1: for K >= 3 and m coprime to K, reduced to 1 <= m <= (K-1)/2, a tuning in
Omega_cv(Q) lies in Omega_{K,m} iff

    A_{K,m} - q_m cot(pi/K) >= 0,   M_m = p_m - i q_m,   A = p(1 - p) - q^2.

Their Appendix B excludes all large K at once: membership forces p in [0, 1], hence cos(2 pi m/K) in an explicit
rational interval on which sqrt(1 - u^2) >= sin_min > 0, and with A <= 1/4 it gives
cot(pi/K)^2 <= R := 1/(16 c^2 sin_min^2), c = (1 - beta)/((L - mu) gamma) with gamma = alpha. R and sin_min^2 are
exact rationals;
K_0 = min{K : cot(pi/K)^2 > R} and every mode with K < K_0 are decided by interval arithmetic.
`omega_aug_report` reproduces Table 1 of Fotopoulos et al. (all 13 rows, 9 digits), their interval
[-56921/74395, 168555301/185987500], sin_min^2 and K_0 = 10 at P_1 exactly (tests/test_regions.py).

**Identity with Pasechnyuk-Vilensky's Theta.** With z = e^{i omega}, a(z) = z^2 - (1 + beta - gamma mu) z + beta and
b(z) = z^2 - (1 + beta - gamma L) z + beta,

    (L - mu) gamma M = -a(z)/z,      (L - mu) gamma (1 - M) = b(z)/z,

so -arg M + arg(1 - M) is the angle from -a(z) to b(z), which is Theta(omega) by Pasechnyuk-Vilensky's
definition. Hence Omega_{K,m} = {Theta(2 pi m/K) <= pi/K}: his single-harmonic region of Theorem 7.1(ii) is exactly
Fotopoulos et al.'s Omega_aug, as he states. The identity is checked symbolically (SymPy) and numerically on
26,700 (tuning, K, m) triples.

## 5. Numerical MILP (`hb1d.milp`)

One binary o_ij per pair (o_ij = 1 means x_i <= x_j), a margin t in [0, 1] and big-M switching of the rows
x_j - x_i >= t, g_j - g_i - mu (x_j - x_i) >= t, L (x_j - x_i) - (g_j - g_i) >= t. Maximizing t gives in one
solve: t* > 0 (a strict cycle), t* = 0 (only degenerate cycles), infeasible, or a time limit (no verdict).

Big-M: with x in [0, 1] (x_0 = 0 is the minimum and sum x = 1), |g_k| <= ((1 + beta) + beta + 1)/alpha =
2(1 + beta)/alpha, so |g_j - g_i| <= 4(1 + beta)/alpha and |x_j - x_i| <= 1; every switched expression plus
t <= 1 is bounded by M = 4(1 + beta)/alpha + L + 1 (we use 1.01 M; about 4.0 at Prop. 9.4). The order rows need
M = 2. Tag: **[numerical]**; the MILP is used only as a cross-check (tests/test_exact.py) because the exact
search is faster at the lengths that matter (HiGHS reaches its 60 s limit at K = 13 for Prop. 9.4).

## 6. Phase diagram (`scripts/phase_diagram.py`)

Grid in Pasechnyuk-Vilensky's unrolled coordinates: rho = rotation number of the tangency point on the boundary of
Omega_1 (where Lc(c_*) = Lc'(c_*) = 0, c_* = cos(2 pi rho)), delta = distance in alpha L beyond that point at
fixed beta, normalized by the height delta_top(rho) at which the vertical line enters the K- or (K+1)-tongue.
Heights h = delta/delta_top run over 25 levels below the gap tops (h < 1) and 4 levels above them
(h = 1.03 ... 1.15, cells in the neighbouring tongues); h = 1 itself is skipped because it is the tongue boundary.
Each grid tuning is rounded to a rational and classified (Omega_1, Omega_aug, mixed-condition lengths). For the
mixed condition, `hb1d.mixed.mixed_float_record` solves both LPs at every K <= 30 (the margin sigma of a feasible
weight vector and the margin tau of an infeasibility vector) and returns +1 if sigma > 1e-8, -1 if tau > 1e-8, and
0 (undecided) otherwise or on any solver failure; the record is stored per cell. The exact 1D search runs at the
smallest Omega_aug length and at every mixed-feasible or undecided length K <= 10. Cell labels are [numerical]
(floating-point classification; Omega_aug membership is rechecked rigorously in every cell by
`scripts/export_rigorous.py`); the 1D searches inside each cell are exact.
