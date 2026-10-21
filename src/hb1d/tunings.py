"""Named heavy-ball tunings used throughout, as exact rationals.

Conventions (shared by every module): L = 1, mu = 1/kappa with kappa = L/mu >= 1 (Pasechnyuk-Vilensky's
convention; Fotopoulos et al. write kappa_F = mu/L = 1/kappa), step size alpha = gamma (alpha L = alpha), and

    x_{k+1} = x_k - alpha f'(x_k) + beta (x_k - x_{k-1}),

the same recursion in both papers.
"""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class Tuning:
    name: str
    alpha: F        # alpha L
    beta: F
    kappa: F        # L / mu

    @property
    def mu(self) -> F:
        return 1 / self.kappa


def polyak(kappa) -> Tuning:
    """Polyak's tuning alpha = 4/(sqrt L + sqrt mu)^2, beta = ((1 - sqrt q)/(1 + sqrt q))^2, q = 1/kappa.
    Rational whenever kappa is a rational square (kappa = 100 gives 400/121 and 81/121)."""
    kappa = F(kappa)
    r = F(int(round(float(kappa) ** 0.5)))
    if r * r != kappa:
        raise ValueError("Polyak's tuning is rational only for square kappa")
    s = 1 / r                                    # sqrt(mu / L)
    return Tuning(f"Polyak kappa={kappa}", 4 / (1 + s) ** 2, ((1 - s) / (1 + s)) ** 2, kappa)


PROP13 = Tuning("Pasechnyuk Prop. 13 (gap I_4, two-harmonic 13-cycle in R^5)", F(70471, 31800), F(6169, 9900), F(100))
CONV_GAP = Tuning("Pasechnyuk Prop. conv-gap (gap I_4, converges)", F(111, 50), F(31, 50), F(100))
POLYAK100 = polyak(100)
# Fotopoulos et al. Example 3.2: kappa_F = mu/L = 4/5, beta = 17/20, L gamma = 18/5; in int(Omega_{5,2}) \ Omega_rou.
FHDA_EX32 = Tuning("Fotopoulos et al. Example 3.2 (in Omega_{5,2}, not in Omega_rou)", F(18, 5), F(17, 20), F(5, 4))
# Fotopoulos et al. point P_1 (Section 5): kappa_F = 531/10^4, beta = 4879/10^4, gamma = 32829/12500; 4D cycle with
# modes 2/10 and 3/10, outside Omega_aug.
FHDA_P1 = Tuning("Fotopoulos et al. P_1 (4D cycle, outside Omega_aug)", F(32829, 12500), F(4879, 10000), F(10000, 531))

# Tunings found here inside Pasechnyuk's first-level fan tongues at kappa = 100 (scripts/find_fan_tunings.py):
# strictly mixed-feasible at level k with modes 1/k, 2/k (so a cycle exists in R^4), certified outside Omega_aug.
FAN7 = Tuning("fan tongue 2/7 in gap I_3 (this repo)", F(34, 13), F(312, 665), F(100))
FAN9 = Tuning("fan tongue 2/9 in gap I_4 (this repo)", F(2201, 986), F(199, 320), F(100))

NAMED = {"fan7": FAN7, "fan9": FAN9, "polyak100": POLYAK100, "prop13": PROP13, "convgap": CONV_GAP, "ex32": FHDA_EX32, "p1": FHDA_P1}
