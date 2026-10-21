"""Phase-diagram geometry of Pasechnyuk-Vilensky (arXiv 2610.09912, Sections 'Regions' and 'Levels').

x = alpha mu, y = alpha L (L = 1). The lag-one quadratic

    Lc(c) = 4 beta c^2 + 2 (y - 1 - 2 beta - beta^2 + beta x) c + (2 - 2x + xy + 2 beta^2 - 2 beta y)

decides Theta(omega) > omega/2 <=> Lc(cos omega) > 0. Omega_1 (lag-one convergence) is {Lc > 0 on [-1, 1]};
the K-tongue (= Omega_{K,1} of Fotopoulos et al. = GTD's root-of-unity region) is {Lc(cos(2 pi/K)) <= 0}.

Unrolled coordinates (Pasechnyuk Fig. regions(d), level1): a boundary point of Omega_1 with tangency rotation number
rho = omega_*/2 pi has Lc(c_*) = Lc'(c_*) = 0 at c_* = cos(2 pi rho); the tuning with unrolled coordinates
(rho, delta) is (alpha_b(rho) + delta, beta_b(rho)), i.e. delta is the distance in alpha L beyond the boundary.
"""
import math

from scipy.optimize import brentq

from .tunings import Tuning


def lag_one(t: Tuning, c):
    al, be, mu = float(t.alpha), float(t.beta), float(t.mu)
    x, y = al * mu, al
    return 4 * be * c * c + 2 * (y - 1 - 2 * be - be * be + be * x) * c + (2 - 2 * x + x * y + 2 * be * be - 2 * be * y)


def boundary_point(kappa, rho):
    """(alpha, beta) on the boundary of Omega_1 with tangency at c = cos(2 pi rho), 0 < rho < 1/2 (floats)."""
    mu = 1.0 / float(kappa)
    c = math.cos(2 * math.pi * rho)

    def alpha_of(be):              # from Lc'(c) = 0
        return ((1 + be) ** 2 - 4 * be * c) / (1 + be * mu)

    def F(be):                     # Lc(c) at Lc'(c) = 0:  -4 beta c^2 + constant term
        al = alpha_of(be)
        x, y = al * mu, al
        return -4 * be * c * c + 2 - 2 * x + x * y + 2 * be * be - 2 * be * y

    grid = [i / 2000 for i in range(1, 2000)]
    for b0, b1 in zip(grid, grid[1:]):
        if F(b0) * F(b1) < 0:
            be = brentq(F, b0, b1, xtol=1e-15)
            return alpha_of(be), be
    raise ValueError(f"no boundary point at rho = {rho}")


def unrolled_to_tuning(kappa, rho, delta):
    al, be = boundary_point(kappa, rho)
    return al + delta, be


def tuning_to_unrolled(t: Tuning):
    """Inverse map at fixed beta: the rho with beta_b(rho) = beta (on the branch 0 < rho < 1/2) and delta."""
    be = float(t.beta)
    f = lambda r: boundary_point(t.kappa, r)[1] - be
    rs = [i / 1000 for i in range(20, 500)]
    for r0, r1 in zip(rs, rs[1:]):
        try:
            if f(r0) * f(r1) < 0:
                r = brentq(f, r0, r1, xtol=1e-14)
                return r, float(t.alpha) - boundary_point(t.kappa, r)[0]
        except ValueError:
            continue
    raise ValueError("not found")


def in_omega1(t: Tuning, n=4001):
    return all(lag_one(t, -1 + 2 * i / (n - 1)) > 0 for i in range(n))


def dip(t: Tuning, n=20001):
    """Rotation-number interval {omega/2pi : Lc(cos omega) < 0} (floating point), or None."""
    pts = [i / (2 * (n - 1)) for i in range(1, n)]
    neg = [r for r in pts if lag_one(t, math.cos(2 * math.pi * r)) < 0]
    return (min(neg), max(neg)) if neg else None
