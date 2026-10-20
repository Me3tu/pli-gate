# C0 rate ladder (written 2026-10-08 before any C0 ladder extraction)

Added in response to review round 5 (item 4): repeat the registered rate ladder
(2026-08-26_rate_scan/REGISTERED_LADDER.md) on the corrected-yield zero-delay
eight-tube stream C0 (/data/dhritman/feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5),
so that the rate scaling is shown on the corrected yield. This is NOT registered in the
sense of the 2026-08 protocol: it is a post hoc replication and is reported as such.

Method: identical to 2026-09-30_referee2/ladder_joint.py + ladder_boot.py
(per-shard epoch compression t_new = t0/k + (t - t0) per eventID group, the published
clinical extractor sim1_extract.clinical_per_shard imported unchanged, per-shard
rng [90900, shard] shared by every rung and rule, twelve rungs k = 0.15 ... 3.00,
joint shard bootstrap B = 4000 with the same resample for all rungs, GLS with the
bootstrap covariance). mu per rung is measured on C0 itself (all shards, and the
first six shards as in mu_ladder.json).

Primary tube set: the six tubes common to all rungs on the discovery stream
(ladder12.json common_tubes), so that the exponents are comparable; all eight tubes
reported as well if every rung gives a positive mean |DiD|.

Expectations stated before running:
E1  mu(k)/k constant within 1 %, and mu(k=1) on C0 equals the discovery value within 2 %
    (same activity and geometry; yield changes the signal, not the prompt rate).
E2  k = 1 rung mean |DiD| over the six tubes agrees with the C0 production value within
    its bootstrap uncertainty.
E3  exponent within the discovery interval [0.973, 1.017] or within 0.05 of 1.
E4  prefactor (mean |DiD| / mu) smaller than on the discovery stream by roughly the
    corrected/discovery artifact ratio (0.79), since the stronger signal lowers the
    randoms-to-signal ratio.
If E3 fails, the paper reports the C0 exponent as measured and drops any claim of
linearity that rests on the discovery stream alone.
