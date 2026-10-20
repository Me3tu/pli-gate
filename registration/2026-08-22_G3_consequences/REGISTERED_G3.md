# G3 consequences — registered before any of the three are computed

The paired bootstrap measured correlation +0.967 to +0.981 between the two
matching rules, so sigma on their difference is 14-19 % of the independent-arm
value. Three things follow. All three use FROZEN data; no new simulation.

## (a) NEMA significance

Frozen: closest 1.6490 +/- 0.0104 ns, all-prompts 1.7194 +/- 0.0111 ns, shift
+0.0703 ns. Independent-arm sigma = hypot(0.0104, 0.0111) = 0.01521 ns, giving
the published **4.6 sigma**.

**Prediction.** The paired sigma is measured, not assumed, by dumping NEMA
per-shard histograms and resampling shards coherently. If NEMA behaves like the
water and disk phantoms the ratio lands near 0.14-0.19 and the significance
grows by roughly its inverse, to **24-33 sigma**. A ratio far outside that band
is itself the finding and is reported as such.

## (b) Pedestal-law statistics with paired errors

V4's identity-line statistics weighted each region by an arm-difference
uncertainty that was 5-7x too large. Redone with paired errors.

**Registered prediction, taken from the decision and not softened if it holds:**
the slope stays near 0.99, but the residual scatter (RMSE 0.70 %) becomes
statistically significant against the smaller errors. That is, chi2/dof rises
well above 1. If that is what comes out, the claim changes from "agreement
within uncertainties" to "the law captures the dominant term and a resolved
residual structure remains", with per-region pulls shown.

Failure modes, both reportable: chi2/dof stays near 1 (the law is complete at
the new precision), or the slope moves off 0.99 (the reweighting changes the
central result, not only its error).

## (c) Does pairing make the D1 estimator swap powered?

V3b compared the windowed-and-corrected estimator against the raw censored-mean
estimator on identical events and got -15.11 +/- 42.88, +5.92 +/- 18.23 and
+28.38 +/- 10.93 ps across three regions. Those errors came from bootstrapping
events independently per estimator. The two estimators see the SAME events, so
the same pairing argument applies.

**Prediction.** Pairing the replicate draw across the two estimators shrinks
sigma substantially. Whether it shrinks by the same 0.14-0.19 factor is NOT
predicted: the correlation between two estimators on one event set is a
different quantity from the correlation between two matching rules, and could
plausibly be higher. If the paired sigma falls below about 4 ps, the 9-16 ps
chain offset becomes testable on existing data and D1 stops being underpowered.

**No tuning.** Nothing is adjusted to make any of the three land.
