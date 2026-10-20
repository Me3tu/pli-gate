# Lineage of the eight-tube numbers

This table was Table S16 of a draft of the supplement. The final supplement
(Notes S1 to S8) moved it here. Note, table, figure and equation numbers below
refer to the final paper and supplement.

The eight-tube phantom enters the analysis through three acquisitions.

- The earlier-yield stream (called the discovery stream in the project records
  and in `registration/`) predates the correction of the I3 convention and
  carries about three quarters of the corrected o-Ps yield (Supplementary
  Note S8). Its mechanism results (timing-resolution invariance, pedestal,
  registered rate replay) are replicated on the corrected-yield arm, and its
  magnitudes are superseded.
- The corrected-yield arm carries every eight-tube magnitude of the main text
  at zero delay and its mechanism tests other than the per-triple truth
  instrumentation.
- The corrected-yield arm and the realistic-timing arm share the macro, binary,
  material file, seeds, shard layout and converter windows and differ only in
  the direct delay.

A dash marks a quantity not evaluated on that dataset. pp: percentage points of
the configured lifetime.

| | earlier-yield stream, 0 ps | corrected yield, 0 ps | realistic timing, 388 ps |
|---|---|---|---|
| build | May 2026 source state, binary not retained (a) | production, `e37be8b4...` | production, `e37be8b4...` (b) |
| role in the paper | replication | zero-delay magnitudes and mechanism tests | realistic-timing magnitudes |
| closest-prompt readout, tubes low by | 15 to 28% (first extraction) | 12.5 to 23.1% | - |
| mean abs. bias, closest -> all-prompts, first extraction | 20.7 -> 6.1% | - | - |
| mean abs. bias, closest -> all-prompts -> a priori acceptance (359 ps), production extractor | 20.5 -> 6.0 -> 1.3% (signed +0.2%) | 18.0 -> 6.7 -> 1.9% (signed -0.7%, Fig. 3b) | - |
| all-prompts with the randoms-arm acceptance, mean abs. bias (signed) | 2.2% (-1.7%) | 2.9% (-2.4%) | - |
| randoms-arm acceptance edge, fit over all regions: width, offset beyond 5 ns | 339 ps, +71 ps | 339 +- 3 ps, +64 ps | - |
| true-triple acceptance edge (Supplementary Note S6) | - | 358.6 +- 2.4 ps, +9 +- 3 ps | - |
| over-subtraction under all-prompts matching, fraction of the true tail counts | consistent with zero | 0.8% (5 sigma, Note S6), on average with the mean delay of the genuine triples within the fit window and no measured lifetime shift (Note S7) | - |
| timing-resolution sweep, tube-mean closest-prompt bias | -17.8 to -17.3% (303 to 50 ps, extractor of that stage) | -18.0 to -17.6% (303 to 50 ps and noiseless, production extractor) | - |
| pedestal plateau, counts per 50 ps bin per tube | 48.9 (20.4 to 22.3% of the randoms estimate) | 49.2 (20.5 to 22.6%) | - |
| template likelihood at 303 ps, mean abs. bias (truncated mean -> template) | 17.8 -> 13.0% | - | - |
| selection shift, inverse-variance mean over the tubes (pp) | 12.22 +- 0.20 | 9.53 +- 0.14 | 4.77 +- 0.05 |
| pedestal law, Eq. (5) | signed regression with disks, water and NEMA: slope 1.032, scatter 0.29 pp | signed regression (Fig. 4c): slope 1.059, scatter 0.39 pp. Out of sample, all eight tubes within 2.0 sigma (0.37 pp scatter), 2.3 sigma with 0.29 pp | - |
| rate dependence | registered replay, exponent 0.997 (0.973 to 1.017), 4.2 pp per 0.01 in mu | replay (c), exponent 0.984 (0.966 to 1.002), 3.3 pp per 0.01 in mu, 0.75 to 0.82 of the earlier-yield value at every rung | - |
| retention of the selection shift, 388 ps over 0 ps | 0.359 +- 0.007 (registered scoring) | 0.437 +- 0.008 | - |
| three-component selection shift, all-prompts pooled bias | - | - | 7.4 to 18.8 pp, -1.0 +- 0.9% (joint likelihood, bootstrap only) |
| event-identity lifetime, corrected minus earlier yield | all events -42 to +20 ps per tube, annihilation within 75 mm of the decay -5.5 to +12.0 ps (chi2 = 3.2 on 8 df, Note S8) | (same comparison) | - |
| event-identity lifetime against configured, annihilation within 75 mm | -0.54 to +0.47% | -0.43 to +0.46% | - |
| known-content test, randoms-correction error under all-prompts matching | - | -1.4 to +0.9% over twelve rates (Note S7) | - |

(a) Inferred from the source timeline. The loader of that state lacks the I3
conversion of Supplementary Note S1.

(b) The in-flight sensitivity arm of Supplementary Note S8 repeats this
acquisition with the variant binary (`cc552a68...`).

(c) Not registered. It used the registered method of the earlier-yield replay
(per-shard epoch compression, published extractor, joint shard bootstrap) on
the six common tubes. Expectations written before the run (rate scaling of mu
within 1%, the k = 1 rung reproducing the production value, exponent within
0.05 of 1, prefactor reduced by about the ratio of the corrected to the
earlier-yield bias) were all met. The record is
`registration/2026-10-08_c0_ladder/REGISTERED_C0_LADDER.md`.

The binary and switches of every run are in `simulation/README.md` (run table)
and `ENVIRONMENT_MANIFEST.json` (`binary_per_run`).
