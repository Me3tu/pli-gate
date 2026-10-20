# Validation

What was checked, with numbers, and what was not. Every figure here is a
measurement from a validated run or a check the test suite reproduces. Nothing
is a target.

Values marked **(reproduced by the suite)** are re-derived every time
`tests/run_tests.sh` runs, so they cannot drift from this document.

## 1. Three-photon final state **(reproduced by the suite)**

From `custom_physics/tests/test_orepowell.cc`, one million sampled events:

    mean photon energy               340.665940 keV   (analytic 340.665940)
    max relative |sum E - 1022 keV|  2.173e-16
    max |sum p| / (m_e c)            1.451e-13
    max ||d| - 1|                    3.331e-16
    photons inside 0 < x <= 1        3000000 / 3000000
    rejection-envelope violations    0

Energy and momentum conserve to machine precision, every photon lands in the
physical region, and the rejection envelope is never violated, so the sampled
spectrum is exact rather than approximate.

Branching, `P(3g) = tau_eff / 142.05 ns`, at six lifetimes:

| tau_eff (ns) | expected | measured | pull |
|---|---|---|---|
| 1.50 | 0.010560 | 0.010530 | −0.42 |
| 1.85 | 0.013024 | 0.013064 | +0.50 |
| 2.00 | 0.014080 | 0.013976 | −1.24 |
| 2.55 | 0.017951 | 0.018040 | +0.94 |
| 2.72 | 0.019148 | 0.018960 | −1.94 |
| 3.39 | 0.023865 | 0.023654 | −1.95 |

All six within 2 sigma, no trend in the pulls.

## 2. Formation channels **(reproduced by the suite)**

From the quickstart run, 200,703 annihilations in water with `I3` = 0.27 in the
o-Ps-intensity convention:

| channel | count | share | expected |
|---|---|---|---|
| direct, 2 gamma | 128,498 | 0.6402 | 0.640 |
| p-Ps, 2 gamma | 18,140 | 0.0904 | 0.090 (= I3/3) |
| o-Ps, 2 gamma | 53,338 | — | — |
| o-Ps, 3 gamma | 727 | — | — |
| **o-Ps total** | 54,065 | **0.2694** | **0.270** (= I3) |

The three-photon share of o-Ps is 727/54,065 = 1.345 %, against
tau/142.05 ns = 1.353 %.

## 3. Lifetime recovery **(reproduced by the suite)**

The quickstart returns **tau_oPs = 1.9072 ± 0.0107 ns** against a configured
1.9215 ns, a difference of −14.3 ps, or −1.3 sigma.

At much higher statistics the sampler was checked against its own input on
5,563,082 annihilations above the 1 ns censor:

    mean(tau | tau > 1 ns) - 1 ns  =  1.83990 ns
    configured                     =  1.84000 ns
    difference                     =  -0.10 +/- 0.78 ps

The physics emits the lifetime it is given, to sub-picosecond.

## 4. Null control

In a single-material phantom, taking the lifetime from the material where the
positron annihilates and from the material where it was created must agree.
Measured difference **+0.00132 ns, +0.07 %** — consistent with zero. In
two-material phantoms deliberately built so the two differ, the same comparison
gives **−3.11 %** and **−2.51 %**, showing the test has the sensitivity to
detect a difference when one exists.

## 5. Over-subtraction (pedestal) law **(reproduced by the suite)**

If a randoms estimator removes more counts than it should, the residue behaves
as a uniform pedestal and biases the recovered lifetime low. Predicting the
difference in bias between two matching rules from the measured difference in
over-subtraction, across 20 regions in 4 phantoms:

    weighted slope against the identity   0.990
    clustered 95 % CI                     [0.948, 1.067]
    RMSE about the identity line          0.70 %

Nothing is fitted; the pedestal mean is the analytic `(T - tau_min)/2`.

The same law was then tested on a configuration it had never seen, with a
finite direct-annihilation lifetime, and predicted it at slope **0.945** with a
*tighter* spread than the original. It is not a property of the prompt-direct
simplification.

## 6. External comparison

Simulated disk and water-cylinder geometries reproduced from published
spectroscopic references were compared against an independent scanner
measurement of the same samples. Agreement is at the **few-percent** level.

**Caveat, stated because it changes how that number should be read**: the
comparison is estimator-dependent. A windowed-and-corrected estimator and a raw
censored-mean estimator applied to the *same* events differ by −15 to +28 ps
depending on the region, because the truncation inversion amplifies statistical
noise. Two of three regions could not resolve a 10 ps effect at all. A
few-percent agreement is therefore a statement about a specific analysis chain,
not about the simulation alone.

---

# Open items

Documented because a user must be able to see what was *not* settled.

## A recovered-lifetime offset of 9 to 16 ps

Recovered lifetimes sit systematically **above** configured ones once the full
reconstruction chain is applied:

    on axis, combined                +9.25 +/- 2.49 ps   (0.50 % of tau)
    off axis, ten positions pooled  +15.62 +/- 1.66 ps   (0.85 % of tau)

The sampler is **not** the cause: it is exact to −0.10 ± 0.78 ps overall, and
remains exact when the analysis' own region masks are applied to its output
(−0.56 ± 4.18, +0.17 ± 1.80, −0.75 ± 1.08 ps in three regions). So the offset
enters downstream of the physics.

An additive per-event delay is **excluded** as the cause, analytically and
numerically. For an exponential censored at c,

    E[t + d | t + d > c] = E[t | t > c - d] + d = c + tau

so the estimator returns tau for any delay d. Numerically, a 50 ps delay moves
the answer by 0.12 ps and a 90.9 ps symmetric timing smear by 0.13 ps. Positron
transport time, photon time of flight and prompt-timing offsets are all
additive and are therefore all excluded. A positive control confirms the test
can see effects it should: a multiplicative scaling and a two-component mixture
are both recovered to sub-picosecond.

**The cause is not known.** It must be something that does not commute with the
censor — a lifetime-dependent selection, a multiplicative scaling, or a change
in population composition.

## An unexplained trend in the calibration closure

Across twelve calibration tubes spanning 1.77 to 2.72 ns, the residual falls
monotonically with configured lifetime:

    slope   -11.70 +/- 5.97 ps/ns      (2.0 sigma from flat)
            -0.448 +/- 0.259  %/ns
    sign change near tau = 1.90 ns

This is **neither** a constant delay (which gives a flat residual in ps) **nor**
a positive fractional bias (which gives one rising with tau). It falls. The
mean signed residual over the twelve is −2.85 ± 2.01 ps, consistent with zero,
so a summary quoted as an RMS or a mean absolute residual does not show this
structure. **The origin is not known**; the additive-delay and timing-smear
causes above are excluded for it too. A constant systematic floor of 3.9 ps,
carrying 23 % of the residual variance, absorbs the remaining scatter.

*Changed in v1.0.1.* Releases up to v1.0.0 quoted −19.11 ± 6.27 ps/ns
(3 sigma), −0.745 ± 0.273 %/ns, a sign change near 2.01 ns and a mean signed
residual of −2.54 ± 2.28 ps. Those came from a superseded pre-correction
analysis pipeline. Every number in this section is now from the
post-correction rerun. No code, test or database changed.

## The early-time model matters more than it looks

Direct annihilation is prompt by default (see `docs/PHYSICS.md`). Setting
`PLI_DIRECT_TAU_PS` to a literature-realistic few hundred picoseconds changes
the answer substantially, because a fraction of direct annihilations then
survives the 1 ns analysis censor:

* the truth-tagged lifetime readout moves from 1.839 to 1.482 ns, −19 %,
  with about 24 % of the fit-window population becoming leaked direct
  annihilation;
* that shift is fully accounted for with no free parameters by the
  count-weighted mean of the component lifetimes, agreeing to 0.07–0.34 %;
* a **windowed-and-truncation-corrected** estimator fares worse, not better,
  on that population: it misses by **−243.8 to −249.8 ps** across regions,
  because the correction assumes a single exponential.

A single-exponential lifetime readout is therefore only meaningful when direct
annihilation really is prompt, or when a mixture model is used instead. This is
the most consequential limitation in this package.
