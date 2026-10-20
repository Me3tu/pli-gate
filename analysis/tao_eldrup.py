"""Tao-Eldrup positronium lifetime model with O2 quenching.

Implements the free-volume pick-off rate from Tao 1972 / Eldrup 1981 and
extends it with the Shibuya 2020 O2 quenching term for PLI simulations
on the Biograph Vision Quadra LAFOV PET/CT.

Unit convention throughout:
  rates in ns^-1, times in ns, radii in nm, O2 partial pressure in mmHg.

This module is the single source of truth for the TE math. It is imported
by:
  - tests/run_tests.sh                  (unit tests)
  - the materials-database generator    (parameter file population)
  - custom_physics JSON consumer        (indirectly, JSON is authoritative
                                         at the C++ side)
"""
import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

LAMBDA_S_NS: float = 2.0
"""Spin-averaged Ps self-annihilation rate (ns^-1). Standard TE value."""

DELTA_R_NM: float = 0.1656
"""Empirical electron-layer thickness (nm). Tao 1972 / Eldrup 1981."""

LAMBDA_3G_NS: float = 1.0 / 142.0
"""Vacuum o-Ps 3-gamma rate (ns^-1) = 1/142 ns ~= 0.00704 ns^-1."""

K_O2_NS_MMHG: float = 0.0204e-3 * 1.26
"""O2 quenching rate constant (ns^-1 mmHg^-1), database v8.

Derived from Shibuya 2020 (Commun Phys 3:173) Eq. 7:
    kappa_O2 = 0.0204 +/- 0.0008 us^-1 umol^-1 L   (at 19 C, per dissolved [O2])
Converted to partial-pressure basis at 37 C with the O2 solubility of human plasma,
    1.26 umol L^-1 mmHg^-1 (Christmas & Bassingthwaighte 2017):
    k_eff = 0.0204e-3 ns^-1 umol^-1 L * 1.26 umol L^-1 mmHg^-1
          = 2.5704e-5 ns^-1 mmHg^-1
Water at 37 C (1.40, Battino 1983) gives 2.856e-5 (sensitivity value). Database v7 and release
v1.0.1 used 2.45e-5, from an unsourced solubility of 1.2 umol L^-1 mmHg^-1 (legacy_v1.0.1/tao_eldrup.py).

docs/PHYSICS.md carries the full derivation chain.
"""


# ---------------------------------------------------------------------------
# Core TE math
# ---------------------------------------------------------------------------

def lambda_po(R_nm: float,
              lambda_S: float = LAMBDA_S_NS,
              delta_R: float = DELTA_R_NM) -> float:
    """Tao-Eldrup pick-off annihilation rate for a spherical free volume.

    ::

        lambda_po(R) = lambda_S * [ 1 - R/R_inf + sin(2 pi R/R_inf) / (2 pi) ]

    where ``R_inf = R + delta_R``.

    Monotonically decreases with R. At R = 0 it equals ``lambda_S`` (Ps is
    fully quenched by pick-off). As R -> infinity it tends to 0.

    Args:
        R_nm:   free-volume radius, nm. Must be non-negative.
        lambda_S: spin-averaged Ps self-annihilation rate, ns^-1.
        delta_R:  electron-layer thickness, nm.

    Returns:
        Pick-off rate, ns^-1.

    Raises:
        ValueError: if R_nm < 0.
    """
    if R_nm < 0:
        raise ValueError(f"R_nm must be non-negative, got {R_nm}")
    R_inf = R_nm + delta_R
    x = R_nm / R_inf  # in [0, 1)
    return lambda_S * (1.0 - x + math.sin(2.0 * math.pi * x) / (2.0 * math.pi))


def tau_eff(R_nm: float,
            pO2_mmHg: float,
            k_O2: float = K_O2_NS_MMHG,
            lambda_3g: float = LAMBDA_3G_NS,
            lambda_S: float = LAMBDA_S_NS,
            delta_R: float = DELTA_R_NM) -> float:
    """Effective o-Ps lifetime (ns) including pick-off + O2 quenching + 3-gamma.

    ::

        tau_eff(R, pO2) = 1 / ( lambda_po(R) + k_O2 * pO2 + lambda_3g )

    Args:
        R_nm:       free-volume radius, nm.
        pO2_mmHg:   oxygen partial pressure, mmHg. Non-negative.
        k_O2:       O2 quenching rate constant, ns^-1 mmHg^-1.
        lambda_3g:  vacuum 3-gamma rate, ns^-1.
        lambda_S:   spin-averaged Ps self-annihilation rate, ns^-1.
        delta_R:    electron-layer thickness, nm.

    Returns:
        Effective o-Ps lifetime, ns.

    Raises:
        ValueError: if pO2_mmHg < 0.
    """
    if pO2_mmHg < 0:
        raise ValueError(f"pO2_mmHg must be non-negative, got {pO2_mmHg}")
    lam = lambda_po(R_nm, lambda_S=lambda_S, delta_R=delta_R)
    lam += k_O2 * pO2_mmHg
    lam += lambda_3g
    return 1.0 / lam


def R_from_tau(tau_ns: float,
               pO2_mmHg: float = 0.0,
               k_O2: float = K_O2_NS_MMHG,
               lambda_3g: float = LAMBDA_3G_NS,
               lambda_S: float = LAMBDA_S_NS,
               delta_R: float = DELTA_R_NM,
               R_lo_nm: float = 1e-6,
               R_hi_nm: float = 2.0,
               xtol: float = 1e-12) -> float:
    """Invert the TE + quenching model to find R given a target lifetime.

    Brent root-find on ``f(R) = tau_eff(R, pO2) - tau_ns``. ``tau_eff`` is
    monotonically increasing in R, so the solution is unique in (0, inf).

    This is how "calibration" materials (biological tissues) get their R
    values: R is back-derived from tau_pub at the assumed pO2. Because of
    that inversion, tissues cannot appear in any "TE predicts tau"
    validation plot -- the model would be predicting its own input.

    Args:
        tau_ns:     target effective lifetime, ns.
        pO2_mmHg:   oxygen partial pressure, mmHg.
        k_O2, lambda_3g, lambda_S, delta_R: see tau_eff().
        R_lo_nm, R_hi_nm: bracket for the root find, nm.
        xtol:       absolute tolerance on R, nm.

    Returns:
        Back-calculated free-volume radius, nm.

    Raises:
        ValueError: if tau_ns is outside the physically achievable range
            for the given pO2 (either below the R=0 lower bound or above
            the R=R_hi_nm upper bound).
    """
    tau_min = 1.0 / (lambda_S + k_O2 * pO2_mmHg + lambda_3g)
    if tau_ns <= tau_min:
        raise ValueError(
            f"tau_ns={tau_ns} ns is at/below the physical minimum "
            f"{tau_min:.6f} ns for pO2={pO2_mmHg} mmHg (would require R<=0)."
        )
    tau_at_hi = tau_eff(R_hi_nm, pO2_mmHg,
                        k_O2=k_O2, lambda_3g=lambda_3g,
                        lambda_S=lambda_S, delta_R=delta_R)
    if tau_ns >= tau_at_hi:
        raise ValueError(
            f"tau_ns={tau_ns} ns exceeds bracket upper limit "
            f"{tau_at_hi:.3f} ns at R={R_hi_nm} nm; widen R_hi_nm."
        )

    def f(R):
        return tau_eff(R, pO2_mmHg,
                       k_O2=k_O2, lambda_3g=lambda_3g,
                       lambda_S=lambda_S, delta_R=delta_R) - tau_ns

    # Bisection, not scipy.optimize.brentq: this module is standard-library only,
    # as README.md states. tau_eff is strictly increasing in R over the bracket,
    # so bisection is safe and the bracket is already validated above. The same
    # pattern is used by analysis/pedestal.py:invert_trunc.
    a, b = R_lo_nm, R_hi_nm
    fa = f(a)
    for _ in range(200):
        mid = 0.5 * (a + b)
        fm = f(mid)
        if fm == 0.0 or 0.5 * (b - a) < xtol:
            return mid
        if (fm < 0.0) == (fa < 0.0):
            a, fa = mid, fm
        else:
            b = mid
    return 0.5 * (a + b)
