# Registered BEFORE any rescaled extraction — the mu ladder

2026-08-26. Zero compute: one dataset already on disk, re-extracted at several
synthetic rates.

## The transformation, and the trap in the obvious version

Rate is density in time, so compressing time raises mu. But a NAIVE global
rescale t -> t/k is WRONG: it compresses the intra-decay delays too, so the
o-Ps lifetime itself would shrink by k and the artefact-vs-mu curve would be
confounded with a tau change. That would test nothing.

The correct operation rescales the DECAY EPOCHS and leaves the within-decay
structure untouched. For every (shard, eventID) group, with t0 the earliest
single in that group:

    t_new = t0 / k  +  (t - t0)

This preserves exactly:
  * every intra-decay delay, so tau is unchanged
  * the number of singles, their energies, positions and truth labels
and changes only the spacing BETWEEN decays, so the accidental-coincidence
rate scales by k and mu scales by k.

Windows are UNCHANGED throughout: 2 ns pair, +/-5 ns prompt, +50 ns randoms.

## The ladder, fixed now

    k     = 0.25, 0.50, 1.00, 1.50, 2.00, 3.00
    mu    ~ 0.0073, 0.0146, 0.0292, 0.0438, 0.0584, 0.0876   (frozen arm, mu_1 = 0.0292)

A 12x span in mu on ONE dataset, same geometry, same seed, same events, fully
paired. Compare with the 2.7x span across four mixed phantoms that gave
n = 1.42.

## Mandatory control

The k = 1.00 rung must reproduce the untouched dataset EXACTLY -- byte-identical
per_voi.csv under both matching rules. If it does not, the rescaling tool is
wrong and nothing on the ladder is trusted. This is the same discipline that
caught the two wrong definitions and that the P5 control failed to provide.

## Predictions

P-L1  k = 1.00 reproduces the untouched extraction byte-for-byte.
P-L2  measured mu scales linearly with k, within 5 %.
P-L3  the artefact scales as mu^n with n determined from the ladder fit;
      registered expectation n = 1.0 to 1.5, i.e. NEARER LINEAR THAN
      QUADRATIC.

P-L3's expectation is stated as the decision's corrected mechanism, not the
original mu^2: the delayed arm removes the first-order accidental CONTENT, and
what survives closest-prompt is the SELECTION bias among kept accidentals,
which is O(mu). The four mixed-phantom points already measured 1.42 / 1.21,
consistent with that. This is a genuine prediction only because the ladder is
a different and much better-controlled measurement; if it returns n near 2 the
corrected mechanism is wrong too.

## What the ladder cannot do

It varies rate at FIXED true activity by editing timestamps. It therefore tests
the analysis-level rate dependence, not detector effects that would accompany a
real activity change (dead time, pile-up), which this simulation does not model
anyway. Stated so the paper does not overclaim it as an activity scan.

## Extension to twelve rungs, registered 2026-08-26 before extraction

Request A1 asks for "the twelve rescaling rungs". Six existed (k = 0.25, 0.50,
1.00, 1.50, 2.00, 3.00). Six are added, chosen to fill the gaps and to probe
the low-mu end where the rules approach degeneracy:

    k added = 0.15, 0.35, 0.70, 1.25, 1.75, 2.50
    full set  0.15 0.25 0.35 0.50 0.70 1.00 1.25 1.50 1.75 2.00 2.50 3.00
    mu span   0.0044 to 0.0877, a factor of 20

Method, controls and windows are unchanged. Predictions:

P-L4  the six added rungs fall on the SAME line as the original six; the
      exponent from all twelve agrees with n = 1.009 within +/- 0.05.
P-L5  k = 0.15 (mu = 0.0044) is at or near the degenerate limit -- the
      diluted C2 run was exactly zero at mu = 0.0040 -- so this rung may
      return a DiD at or below numerical resolution. If it does it is
      EXCLUDED from the fit as a bound, not forced onto the log axis, and
      that exclusion is reported.
