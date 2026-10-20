# Registered before the three-component fits ran — 2026-08-26

## What is being measured

The closest-vs-all-prompts rule shift on the eight Sim-1 tubes under the
THREE-COMPONENT + acceptance estimator, on the existing C2 (388 ps)
histograms. Fitting only; no simulation.

This exists because the truncated mean is not the clinically interpreted
estimator (reviewer item 47). A2 currently shows the mechanism persists under
the truncated mean; this shows whether the quantitative claim survives the
proper estimator.

## Two things fixed before any fit

1. **The acceptance width is Sim-1's own.** sigma_e = 359.0 ps, derived a
   priori from the source geometry in `sigma_geo_sim1.json`. C1 used 312.4 ps,
   which belongs to the Mercolli/water cylinder. Reusing it here would repeat
   the exact error the acceptance memo records.

   A registration already scored in that file: I predicted Sim-1's sigma_geo
   would be SMALLER than NEMA's 351.4 ps because the tubes are axially short.
   **REFUTED** -- it is 353.2 ps, marginally LARGER. The radial offset of the
   80 mm ring matters more than the axial extent. Recorded as a failed
   prediction; the measured value is used regardless.

2. **The estimator is unchanged otherwise.** TAU_PPS = 0.125 ns and
   TAU_DIR = 0.388 ns stay fixed at their b2_estimator values, the Dirichlet
   prior is unchanged, and the fit is the same four-start Nelder-Mead.
   Because TAU_DIR is pinned at 388 ps this estimator is defined for the
   388 ps arm ONLY. The 0 ps arm is NOT fitted here: changing TAU_DIR would
   make it a different estimator and the comparison would not be like-for-like.

## Registered expectation, from the authors' decision

> a rule shift of order 2-4 pp that survives the proper estimator

stated loosely and in advance. NEMA's own component-estimator shift is
1.7123 -> 1.7624 ns, i.e. 0.0501 ns = 2.8 % of the configured 1.80 ns.

    P-C2-1  the mean rule shift over the eight tubes is in 2-4 pp of
            configured tau, and is POSITIVE (all-prompts above closest-prompt,
            i.e. closest-prompt still reads low).
    P-C2-2  the shift is resolved: mean / paired sigma > 3 on the pooled value.
    P-C2-3  the shift is SMALLER than the truncated-mean DiD at 388 ps
            (mean 4.98 pp), because the component estimator models the direct
            term explicitly instead of absorbing it into a single exponential.

If the shift comes out near zero, the quantitative realistic claim does NOT
survive the proper estimator and A2 must be labelled a mechanism-persistence
figure only. That outcome is reportable and will be reported.

## Errors

PAIRED: one shard resample, BOTH rules refitted from exactly those shards.
600 draws. Point estimates use the full four-start fit; replicates use a warm
start from the point estimate, which reproduces it to 0.00 ps on the test case.
