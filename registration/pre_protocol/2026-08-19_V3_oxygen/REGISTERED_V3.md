# V3 — registered before the materials are built and the run is launched

## What the request asked for

"Simple geometry, fixed R, high statistics per point, pO2 in {0, 2, 6, 10, 20,
40, 100} mmHg. Extract tau per point, overlay the analytic
1/tau = 1/tau_TE + lambda_3g + kappa*pO2 line. No fitting: the analytic line is
the prediction. Deliverable: one figure + per-point residuals. This validates
the oxygen implementation independent of clinical detectability."

## What a GATE run can and cannot test here — stated before running

`PLIMaterialLoader.hh:12-14` says it outright: *"Python is the single source of
truth for the Tao-Eldrup + O2 quenching math. The JSON carries the
already-computed tau_TE_ns per material; this C++ loader does NOT re-derive R
or tau."* Confirmed in the source: `pO2_mmHg` is read at
`PLIMaterialLoader.cc:163` and thereafter appears only in a diagnostic print at
line 284.

So the chain splits, and only one half is a GATE question:

| link | tested by |
|---|---|
| kappa_O2 and the Henry's-law conversion from Shibuya 2020 | V10 unit test, arithmetic, **not GATE** |
| `tau_eff(R, pO2)` -> JSON `tau_TE_ns` | Python, checkable exactly |
| JSON -> C++ loader -> sampler -> emitted delays -> extracted tau | **this run** |

This run therefore validates that the pipeline *delivers the lifetime the oxygen
model computed*, end to end, at ten different values of it. It does not and
cannot validate the oxygen physics itself. Reporting it as "the oxygen
implementation is validated" without that split would overstate it.

## The arithmetic that sets the design, computed before the run

R is fixed at **0.272081 nm**, the radius consistent with the existing
`WaterSteinberger` entry (tau = 1.84 ns at pO2 = 0). Then
lambda_po(R) = 0.53643601 /ns, lambda_3g = 0.00704225 /ns,
kappa = 2.45e-5 /ns/mmHg.

| pO2 (mmHg) | tau (ns) | shift vs pO2=0 |
|---|---|---|
| 0 | 1.840000 | 0 |
| 2 | 1.839834 | -0.17 ps |
| 6 | 1.839502 | -0.50 ps |
| 10 | 1.839171 | -0.83 ps |
| 20 | 1.838343 | -1.66 ps |
| 40 | 1.836688 | -3.31 ps |
| 100 | 1.831743 | -8.26 ps |

**The whole physiological grid spans 8.26 ps, 0.449 % of tau.** With the oracle
estimator's sigma = tau/sqrt(n), resolving that span at 3 sigma needs
**4.5e5 genuine triples per point**; at 1 sigma, 5.0e4. This is the same wall
V6 hit from the other side: V6 measured a slope of +1.4 [-5.7, +8.5] ps/mmHg
against a configured -0.0903 ps/mmHg, i.e. the configured effect sits two
orders of magnitude inside the noise.

**A seven-point physiological sweep is therefore a null test by construction,
not a validation.** Running it and reporting "consistent with the analytic
line" would be true and nearly content-free.

## Design, and the deviation from the request, stated as a deviation

I am adding **three non-physiological points at 300, 1000 and 3000 mmHg**
(predicted tau 1.815448, 1.760631, 1.620803 ns; shifts -24.6, -79.4, -219.2 ps).
These make the implementation's response resolvable at achievable statistics
and turn the run into an actual test of the chain. They are labelled
non-physiological everywhere they appear and are not offered as clinical
conditions. The seven requested points are run unchanged and reported in full.

All ten materials are simulated in **one run**, as ten identical 20 mm x 300 mm
cylinders on a 5 x 2 grid at x in {-160, -80, 0, +80, +160} mm, y = +-60 mm,
each its own VOI. One run means one binary, one seed stream and one geometry
for all ten points, so anything that drifts drifts commonly and cancels in the
slope. It does not reduce CPU — cost scales with total decays either way.

Direct annihilation is at **0 ps**: `PLI_DIRECT_TAU_PS` is unset, so this is
the production model, not V2's arm B.

## Predictions

1. **Every point lies on the analytic line within its own statistical error.**
   There is no fitting anywhere; the line is computed from R, kappa, lambda_3g
   and pO2 alone. A systematic offset common to all ten points would indicate
   a JSON->C++ transport error; a *sloped* residual would indicate that kappa
   is being applied differently at the two ends.
2. **The seven physiological points will not separate from each other.** I
   expect their mutual differences to be inside 1 sigma at the statistics I
   can afford. This is a prediction about the experiment, not about the
   physics, and it is why points 8-10 exist.
3. **The three extended points will separate cleanly**, at high significance,
   and their residuals against the analytic line will be the number that
   carries the validation.

## Failure modes, all reportable

* all ten on the line -> chain validated across a 12 % lifetime range; say so
  with the residual RMS, and say plainly that the O2 physics was not what was
  tested
* common offset -> transport bug between JSON and the sampler; outranks the
  rest of V3
* sloped residual -> the extended points are being handled differently from the
  physiological ones; report and stop, do not adjust kappa

**No tuning.** R, kappa, lambda_3g and lambda_S are the frozen project
constants. If a point misses, it is reported as missing.
