"""Membership in the augmented root-of-unity region Omega_aug of Fotopoulos, Helcig, Deck and Alsbach
(arXiv 2610.25104), and Pasechniuk-Vilensky's angle function Theta (arXiv 2610.09912).

Fotopoulos et al., Lemma 3.1: for K >= 3, m in U_K = {1 <= m < K : gcd(m, K) = 1}, after replacing m by K - m
so that 1 <= m <= (K-1)/2, write M_m = p_m - i q_m with

    p_m = ((1 + beta - mu gamma) - (1 + beta) cos(2 pi m/K)) / ((L - mu) gamma),
    q_m = (1 - beta) sin(2 pi m/K) / ((L - mu) gamma).

Then (gamma, beta) in Omega_{K,m}  <=>  (gamma, beta) in Omega_cv(Q)  and  phi_{K,m}(Delta_m) <= 0, which (by
their eq. (3.5) with m Delta_m = 1 mod K) reads

    A_{K,m} - q_m cot(pi/K) >= 0,        A_{K,m} = p_m (1 - p_m) - q_m^2.

Omega_aug is the union over all K >= 3 and m in U_K. Their Appendix B excludes all large K at once:
membership forces A >= 0, hence p in [0, 1], hence u = cos(2 pi m/K) lies in an explicit interval on which
sqrt(1 - u^2) >= sin_min > 0; with A <= 1/4 this gives cot(pi/K) <= 1/(4 c sin_min), c = (1 - beta)/((L-mu) gamma).

Signs are decided with mpmath interval arithmetic (outward rounding), so every "certified" sign below is rigorous;
inputs are exact rationals.
"""
from contextlib import contextmanager
from fractions import Fraction as F
from math import gcd

from mpmath import iv, mp

from .tunings import Tuning

PREC = 256


@contextmanager
def _prec(bits):
    old = iv.prec
    iv.prec = bits
    try:
        yield
    finally:
        iv.prec = old


def _iv(x):
    x = F(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def sign(interval):
    """+1 / -1 if the interval is certified positive / negative, 0 if it contains 0 (undecided)."""
    if interval.a > 0:
        return 1
    if interval.b < 0:
        return -1
    return 0


def in_cv(t: Tuning) -> bool:
    """Omega_cv(Q_{mu,L}) = {|beta| < 1, 0 < gamma < 2(1 + beta)/L}: HB converges on all quadratics."""
    return abs(t.beta) < 1 and 0 < t.alpha < 2 * (1 + t.beta)


def reduce_mode(K: int, m: int) -> int:
    if not (K >= 3 and 1 <= m < K and gcd(m, K) == 1):
        raise ValueError(f"(K, m) = ({K}, {m}) is not a primitive mode")
    return m if 2 * m <= K - 1 else K - m


def pure_mode(t: Tuning, K: int, m: int, prec: int = PREC):
    """Interval values (A_{K,m}, q_m cot(pi/K), A - q cot(pi/K)) of Fotopoulos et al. Lemma 3.1 / Table 1.
    (gamma, beta) in Omega_{K,m} iff the last value is >= 0 (and the tuning is in Omega_cv)."""
    m = reduce_mode(K, m)
    with _prec(prec):
        a, b, mu = _iv(t.alpha), _iv(t.beta), _iv(t.mu)
        den = (1 - mu) * a
        ang = 2 * iv.pi * m / K
        p = ((1 + b - mu * a) - (1 + b) * iv.cos(ang)) / den
        q = (1 - b) * iv.sin(ang) / den
        A = p * (1 - p) - q * q
        qc = q * iv.cot(iv.pi / K)
        return A, qc, A - qc


def large_K_bound(t: Tuning):
    """Appendix B of Fotopoulos et al.: an exact rational R such that a pure mode at level K can only be feasible
    if cot(pi/K)^2 <= R. Returns (R, u_lo, u_hi, sin_min_sq) with all entries exact Fractions."""
    a, b, mu = t.alpha, t.beta, t.mu
    den = (1 - mu) * a
    P0, P1 = (1 + b - mu * a) / den, (1 + b) / den      # p(u) = P0 - P1 u, P1 > 0
    lo, hi = max(F(-1), (P0 - 1) / P1), min(F(1), P0 / P1)   # p in [0, 1]  <=>  u in [lo, hi]
    if lo > hi:
        return F(0), lo, hi, None                         # p never in [0, 1]: no pure mode at all
    sin_min_sq = 1 - max(lo * lo, hi * hi)
    if sin_min_sq <= 0:
        raise ValueError("u-interval reaches +-1; the uniform bound does not apply (tuning outside Omega_cv?)")
    c = (1 - b) / den
    return 1 / (16 * c * c * sin_min_sq), lo, hi, sin_min_sq


def first_excluded_K(R: F, prec: int = PREC) -> int:
    """Smallest K >= 3 with cot(pi/K)^2 > R, certified; cot(pi/K) increases with K, so all K >= K0 are excluded."""
    K = 3
    with _prec(prec):
        Ri = _iv(R)
        while sign(iv.cot(iv.pi / K) ** 2 - Ri) != 1:
            K += 1
    return K


def omega_aug_report(t: Tuning, prec: int = PREC):
    """Decide membership of a rational tuning in Omega_aug with a rigorous certificate.

    Returns a dict with 'member' in {True, False, None}: True with a witness (K, m) whose value is certified >= 0
    (> 0) or exactly on the boundary; False with K0 and, for every primitive mode with K < K0, a certified
    negative value; None if some sign could not be decided at this precision."""
    out = {"tuning": t.name, "in_cv": in_cv(t)}
    if not out["in_cv"]:
        out["member"] = False
        out["reason"] = "outside Omega_cv(Q): Omega_aug is defined inside Omega_cv"
        return out
    R, lo, hi, smin = large_K_bound(t)
    K0 = first_excluded_K(R, prec)
    out.update(R=R, u_interval=(lo, hi), sin_min_sq=smin, K0=K0)
    modes, undecided, witnesses = [], [], []
    for K in range(3, K0):
        for m in range(1, (K - 1) // 2 + 1):
            if gcd(m, K) != 1:
                continue
            A, qc, v = pure_mode(t, K, m, prec)
            s = sign(v)
            modes.append((K, m, float(mp.mpf(v.mid)), s))
            if s == 1:
                witnesses.append((K, m))
            elif s == 0:
                undecided.append((K, m))
    out["modes"] = modes
    out["witnesses"] = witnesses
    out["undecided"] = undecided
    out["member"] = True if witnesses else (None if undecided else False)
    return out


# ----------------------------------------------------------------------------------------------------------------
# Pasechniuk-Vilensky's angle function, computed from its definition and from the closed form, and the identity
# with Fotopoulos et al.'s angle.

def theta_def(t: Tuning, omega):
    """Theta(omega) in (0, pi]: the angle from -a(e^{i omega}) to b(e^{i omega}) (Pasechniuk Def. 'Angle function'),
    a(z) = z^2 - (1 + beta - alpha mu) z + beta, b(z) = z^2 - (1 + beta - alpha L) z + beta. Floating point."""
    import cmath
    z = cmath.exp(1j * omega)
    a = z * z - (1 + float(t.beta) - float(t.alpha * t.mu)) * z + float(t.beta)
    b = z * z - (1 + float(t.beta) - float(t.alpha)) * z + float(t.beta)
    return cmath.phase(b * (-a).conjugate())


def theta_closed(t: Tuning, omega):
    """Theta from Pasechniuk Lemma 'angle at the threshold': cot Theta = Xi(cos omega)/sin omega."""
    import math
    al, be, mu = float(t.alpha), float(t.beta), float(t.mu)
    A, B = 1 + be - al * (1 + mu) / 2, al * (1 - mu) / 2
    c = math.cos(omega)
    xi = (B * B - (1 - be) ** 2 - A * A + 2 * A * (1 + be) * c - 4 * be * c * c) / (2 * B * (1 - be))
    return math.atan2(math.sin(omega), xi)        # Theta in (0, pi) with cot Theta = xi / sin omega


def fhda_angle(t: Tuning, K: int, m: int):
    """-arg M_m + arg(1 - M_m) of Fotopoulos et al. Lemma 3.1(ii) (m reduced to <= (K-1)/2). Floating point."""
    import cmath
    m = reduce_mode(K, m)
    al, be, mu = float(t.alpha), float(t.beta), float(t.mu)
    w = cmath.exp(2j * cmath.pi * m / K)
    M = ((1 + be - mu * al) - w - be / w) / ((1 - mu) * al)
    return -cmath.phase(M) + cmath.phase(1 - M)
