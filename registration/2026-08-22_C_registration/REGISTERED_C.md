# C-chain — predictions registered BEFORE the realistic-timing runs exist

Written 2026-08-22 ~19:15. The C-chain cannot run in this block: it needs the
cores, and the Mercolli exposure run holds them until ~16:00 on 23 August, about
six hours after this block ends. These predictions are therefore recorded now, so
that when C1 and C2 do run they are tested against a bar written before any of
their data existed.

## The basis

V2 measured the same thing on water: the closest-versus-all-prompts artefact
under direct annihilation at 0 ps, and again at 388 ps, on identical events.
Its registered statistic is the difference-in-differences, and the quantity
carried forward here is the RETAINED FRACTION, that is

    retained = |artefact at 388 ps| / |artefact at 0 ps|

measured on four water regions as 0.561 (bg), 0.737 (VOI_5x15), 0.589 (VOI_10x20)
and 0.576 (VOI_15x25) -- a band of
**0.56 to 0.74**.

## Predictions, per configuration

P-C1  **NEMA IQ under PLI_DIRECT_TAU_PS=388, both rules.** The absolute
      closest-vs-all-prompts artefact retains 0.56 to 0.74 of its 0 ps value.
      The 0 ps value is frozen: closest 1.6490 ns, all-prompts 1.7194 ns, shift
      +0.0703 ns. So the predicted 388 ps shift is **+0.039 to +0.052 ns**.

P-C2  **Sim-1, eight tubes, under 388, both rules.** The same band applies per
      tube. Sim-1 spans configured tau from 1.77 to 2.72 ns, so if the retained
      fraction turns out to depend on tau, this is the configuration that shows
      it. A tau-dependence would be a NEW result, not a failure of the band.

P-C3  **The truncated-mean estimator fails visibly at 388 ps in every
      configuration.** From B5's water measurement the recovered lifetime falls
      by 526 ps under closest-prompt and 564 ps under all-prompts. The
      prediction is a fall of the same order, 0.4 to 0.7 ns, in NEMA and Sim-1
      too. The three-component estimator, which passed its gate today, carries
      the actual number.

## Falsifiers, stated now

- A retained fraction OUTSIDE 0.56-0.74 in a configuration means the 388 ps
  effect is not a simple rescaling of the artefact, and the water measurement
  does not generalise. Report it, do not widen the band.
- A truncated-mean fall much smaller than 0.4 ns would mean water is not
  representative of the leakage into the fit window, which is a statement about
  the phantom rather than about the estimator.
- Every rule difference gets the PAIRED bootstrap. An independent one is the
  error this project has already made once.

---

# Per-configuration numbers, fixed 2026-08-24 BEFORE either C run is launched

The band 0.56-0.74 is the retained ABSOLUTE artefact fraction measured on V2
water. Applied to each configuration's own frozen 0 ps values.

## C1 — NEMA IQ at PLI_DIRECT_TAU_PS=388, both rules

Frozen 0 ps (freeze `nema_closest_vs_allprompts`): closest 1.6490 ns,
all-prompts 1.7194 ns, shift **+0.0703 ns**, paired sigma 2.05 ps.

    predicted 388 ps shift = 0.56 to 0.74 x 0.0703  =  +0.039 to +0.052 ns

NEMA is the sharpest test of the band because it has ZERO lifetime contrast:
one material throughout, so nothing but the rule can move the two arms apart.

## C2 — Sim-1 eight tubes at 388, both rules

Frozen 0 ps DiD per tube (`delta_test.json`), in % of tau:

| tube | tau_cfg | DiD at 0 ps | predicted at 388 ps |
|---|---|---|---|
| T7_brain_tumor | 1.770 | 11.05 | 6.19 - 8.18 |
| T8_myxoma | 1.920 | 15.70 | 8.79 - 11.62 |
| T1_water | 1.921 | 12.13 | 6.79 - 8.98 |
| T5_muscle | 2.030 | 19.58 | 10.96 - 14.49 |
| T4_hepatic | 2.040 | 16.39 | 9.18 - 12.13 |
| T2_adipose | 2.540 | 14.44 | 8.09 - 10.69 |
| T3_adipose | 2.720 | 9.63 | 5.39 - 7.13 |
| T6_brain_healthy | 2.720 | 9.82 | 5.50 - 7.27 |

Per-tube bands, not a pooled one: the prediction is that EACH tube retains
0.56-0.74 of its own 0 ps DiD.

## What each outcome means

- **Inside the band, per configuration.** The 388 ps effect is a simple
  rescaling of the artefact and V2's water measurement generalises.
- **Outside the band.** It does not generalise. Report the measured retention,
  do NOT widen the band retrospectively.
- **A tau-dependence of the retention across Sim-1's eight tubes** (1.77 to
  2.72 ns) would be a NEW result, not a failure. Sim-1 is the only
  configuration that can show it; NEMA cannot, having one material.

## Method, fixed now

Both estimators on IDENTICAL events: the production truncated mean, and the
acceptance-corrected three-component fit that passed its gate. Every rule
difference gets the PAIRED (shard-coherent) bootstrap. Quoted uncertainties
carry the endpoint-sweep systematic in quadrature, as ruled for B4.
