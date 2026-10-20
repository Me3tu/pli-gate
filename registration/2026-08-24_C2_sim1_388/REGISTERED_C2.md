# C2 — registered BEFORE launch, 2026-08-24

Two things are registered here: the original one-factor band (unchanged, from
REGISTERED_C.md), and the decision's post-hoc two-factor model, which earns its
place only by predicting C2 before C2 exists.

## 1. The original one-factor prediction (unchanged, still on the record)

Each tube retains 0.56-0.74 of its own 0 ps DiD. Per-tube predictions are in
REGISTERED_C.md and are NOT restated or adjusted here.

C1 already failed this on NEMA: retention 0.4655 against 0.56-0.74, 9.4 sigma
below the lower edge. That failure stands permanently whatever C2 shows.

## 2. The decision's two-factor model (POST HOC, labelled as such everywhere)

    retention_i = (phi_1/phi_0)_i  x  lever(tau_388,i)/lever(tau_0,i)

with phi measured PER TUBE from C2's own randoms arm (NOT borrowed from
water), and lever computed on the same W = 4 ns window for both arms.

### Predictions

C2-a  Retention falls monotonically with tau_0 across the eight tubes.
C2-b  The eight measured retentions agree with the two-factor model to
      within 15 % RMS.
C2-c  The one-factor (absolute-transfer) model fails on at least half the
      tubes, i.e. >= 4 of 8 fall outside 0.56-0.74.

If C2-b fails, both the two-factor model and the original band are dead and
the retention is configuration-specific with no closed form. That is itself
the reportable result and is accepted in advance as such.

### UNRESOLVED — the lever definition does not reproduce the decision's number

The decision states lever = (m_ped - m1)/m1' and quotes
lever(1.137)/lever(1.684) = 0.806 for NEMA.

The NEMA arm-averages reproduce exactly: 0 ps (1.6490+1.7194)/2 = 1.6842,
388 ps (1.1202+1.1529)/2 = 1.1366, matching the quoted 1.684 -> 1.137. So the
tau inputs are agreed.

The lever ratio does not. With m_ped = 2.0 ns and the censored-exponential
first moment m1 = tau_min + tau (tau_min = 1.0), so m1' = dm1/dtau = 1:

    (m_ped - m1)/m1'      ->  0.1996     (the literal stated formula)
    (m_ped - tau)         ->  2.7342
    (m_ped - tau)/tau     ->  4.0516
    tau                   ->  0.6748
    tau/(m_ped - tau)     ->  0.2468
    sqrt(tau)             ->  0.8215     (closest, but not the stated formula)

None gives 0.806. The exponent that would, if lever = tau^a, is a = 0.5484,
which is not a form anyone wrote down.

**This is not being resolved by choosing whichever function lands nearest
0.806.** That is precisely how the route-1 acceptance error was made, and the
decision on that occasion was that it was a specification error, not a physics
failure. The lever definition needs confirming before C2-b can be scored.

Until it is confirmed, C2-a and C2-c are scoreable and C2-b is held.

## 3. Frozen inputs, fixed now so they cannot drift

Per-tube 0 ps DiD and phi are already frozen in
2026-08-17_derivation_test/delta_test.json and are restated here as the
baseline C2 is measured against:

    tube                 tau_cfg   DiD_0ps %   phi_closest   phi_all
    T7_brain_tumor        1.770      11.05       0.09105      0.00316
    T8_myxoma             1.920      15.70       0.13801     -0.00516
    T1_water              1.921      12.13       0.10574      0.00461
    T5_muscle             2.030      19.58       0.19443      0.01355
    T4_hepatic            2.040      16.39       0.19766     -0.01342
    T2_adipose            2.540      14.44       0.14332      0.00657
    T3_adipose            2.720       9.63       0.08773     -0.00297
    T6_brain_healthy      2.720       9.82       0.09162      0.00424

Note for scoring C2-a: phi_0 is NOT monotone in tau_cfg (it peaks at
T4/T5 around 2.03-2.04 ns and falls at both ends), so a monotone retention
prediction is a genuine constraint on the lever term rather than a restatement
of phi.

## 4. Method

Same macro, activity, digitizer and thread pattern as the frozen 0 ps arm; the
two arms differ ONLY in PLI_DIRECT_TAU_PS. Both matching rules, both
estimators, PAIRED (shard-coherent) bootstrap on every rule difference.

CONVERTER: keep-511 435 870, keep-prompt 870 1300. Union is the full digitizer
band [435, 1300] with no gap and no overlap, so nothing is discarded. Tight
flags unchanged at [481,541] and [1100,1200] so the comparison with the frozen
arm is exact. This is the decision's requirement and the direct consequence of
retraction 2: the ROOT files die at conversion and there is no second chance.

---

# APPENDED 2026-08-25, before the rate-matched re-run produced anything

## Why there is a re-run

The first 128-shard run matched the frozen arm on primaries (250,000,000) and
on singles/primary (1.267 vs 1.272) but NOT on duration. Its shards spanned
1.0000 s against the frozen arm's 0.1382 s, so the instantaneous rate was 7.25x
too low and prompts almost never competed. Measured:

    mu = prompts per +/-5 ns window     low-rate run  0.004
                                        frozen arm    0.029

At mu = 0.004 the closest and all-prompts rules select the same single
candidate, so the DiD was identically 0.00 on all eight tubes. That run is NOT
discarded: it is the low-rate arm of a rate scan (2026-08-25_C2_lowrate_null,
data at /data/dhritman/c2_sim1_388_LOWRATE_null).

## The three-part probe check, now mandatory

A matched run must agree on primaries, duration AND the resulting rate. The
rate is the one the physics depends on. Probe of the corrected configuration,
one shard, run to completion:

    singles per primary   1.2663   (frozen 1.272)
    duration              0.137933 s  (frozen 0.1382)
    mean prompt gap       343.0 ns    (frozen 342.7)
    mu                    0.0292      (frozen 0.029)

## Registered predictions, BEFORE the corrected run

P-R1  The corrected run reproduces mu = 0.029 within +/- 15 %.

P-R2  The artefact returns to the frozen arm's scale: per-tube DiD falls
      inside the per-tube bands registered on 24 August.

P-R3  Across every rate point available, the artefact scales as mu^n with
      n = 2 within +/- 0.5.

Basis for n = 2: the rules differ only when a pair carries more than one prompt
in its window, so the fraction of AFFECTED pairs scales with mu; the
over-subtracted pedestal scales with the accidental rate as well, so phi also
scales with mu. Their product gives mu^2.

Consistency check already available: mu = 0.004 against 0.029 is a factor 53
suppression under mu^2, so the frozen 14.6 pp DiD would appear as 0.28 pp --
consistent with the 0.00 measured. Linear scaling would have predicted 2.0 pp,
which would have been visible. The low-rate null therefore already favours
n = 2 over n = 1, and P-R3 is the falsifiable test of it.

P-R3 is NOT a repair of the failed 24 August band. That band failed because it
transferred a rate-specific number between configurations at different rates;
P-R3 is the quantitative law that explains why, and its exponent is the test.

## Configuration of the corrected run

    setTotalNumberOfPrimaries  1,953,125 per shard   (runner, primaries mode)
    setTimeStop = setTimeSlice 0.137933 s per shard  (= 1,953,125 / 14.16 MBq)
    shards                     128  (10 x 12 + 1 x 8)
    PLI_DIRECT_TAU_PS          388
    converter keep-window      435-870 / 870-1300 keV, unchanged
    tight windows              481-541 / 1100-1200 keV, unchanged
