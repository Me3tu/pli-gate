# R1 — registered before the test

## The question

V3b showed the sampler is exact (−0.10 ± 0.78 ps) and the reconstruction adds
+9 to +16 ps. The twelve-tube closure does not show it (−2.54 ± 2.28 ps). The
two use different estimators:

    oracle (V3/V3b) : mean(tau | tau > 1 ns) - 1 ns, NO upper cut, no correction
    tubes           : same but WINDOWED at W_PROMPT = 5 ns, then truncation_correct

R1 asks whether the offset appears under the tube estimator.

## Design: paired, on identical events

Both estimators are applied to **the same per-event delay array**, so the
comparison between them carries almost no statistical error — they are perfectly
correlated by construction. This is why the test needs no new simulation and no
high-statistics run: it measures a deterministic transformation, not a mean.

Two datasets already on disk, giving two different true lifetimes:

* `2026-06-19_steinberger_validation/water_singles_truetime.h5` — production
  model, direct at 0 ps, true tau ~ 1.84 ns
* `2026-08-19_V2_earlytime/singles_water_388ps.h5` — arm B, whose fit-window
  population has an effective tau ~ 1.48 ns

Two points is not a tau-dependence measurement, but it is two independent checks
of whether the estimator swap moves the answer at all.

## Predictions

1. **The two estimators will NOT agree.** The upper cut at 5 ns discards the
   tail, biasing the windowed mean low; `truncation_correct` exists to undo
   exactly that, so the corrected value should return close to the raw one. The
   question is whether "close" is closer than the 9-16 ps at issue.
2. **If the tube estimator reproduces the raw oracle to within a few ps**, then
   the estimator is NOT what separates the tubes from V3/V3b, the thread splits,
   and the tube closure and the reconstruction offset are two findings.
3. **If the tube estimator differs by roughly the offset**, the thread is one
   finding and D2's −19.11 ± 6.27 ps/ns trend should follow from propagating a
   constant per-event delay through the same inversion. That propagation is the
   second half of R1 and is computed analytically, not fitted.

**No tuning.** The estimator, window, binning and correction are the production
ones, imported rather than re-implemented.
