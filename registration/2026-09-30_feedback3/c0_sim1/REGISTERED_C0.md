# C0: corrected-physics zero-delay eight-tube arm, registered BEFORE the run (2026-09-30 21:50, arch clock; file mtime is the record)

Why. The frozen zero-delay eight-tube stream (2026-06-10_brain_tumor_1B_truetime/sim1_singles_truetime.h5)
predates the I3-convention correction of 2026-08-14. A count test on existing files shows that in that
stream every tube carries 3/4 of the corrected o-Ps yield: the truth-signal count in [1,5] ns of C2
(388 ps, corrected binary) over the frozen stream is 1.74-2.68 per tube, against 1.29-1.93 predicted if
only the direct delay differed and 1.71-2.57 predicted if the frozen stream also had o-Ps yield x3/4
(all eight tubes match the second model within 2-5%). The manuscript states that the 0 ps and 388 ps
arms differ only in the direct delay (Sec. 3.9, Fig. 11, Table S14, Fig. 10b); that is not true for the
eight-tube phantom. NEMA, disks and water were simulated after the correction and are unaffected.

Arm. C2 exactly (2026-08-24_C2_sim1_388): same macro sim1_sc44_388_timedriven.mac, production binary
(software/install/gate-st, sha256 e37be8b4...), materials/pli_materials.json, time-driven,
10 batches x 12 shards + 1 x 8, 0.137933 s per shard, same runner, same converter windows
(keep [435,870]/[870,1300], tight [481,541]/[1100,1200]). Only difference from C2: PLI_DIRECT_TAU_PS
unset (direct annihilation instantaneous, the default). Output /data/dhritman/feedback3_c0/c0_sim1_0ps.

Predictions (numbers from the count test and the pedestal law, written before any C0 data exist):
Q1  truth-signal count in [1,5] ns, C0 / frozen stream: 1.25-1.42 per tube (o-Ps yield x4/3 plus p-Ps).
Q2  truth-grouped (oracle) lifetimes: C0 minus frozen consistent with zero per tube (|z| < 3) and pooled
    (|z| < 2); the lifetime does not depend on the yield.
Q3  prompt competition mu: C0 within 0.5% of C2 (same activity, same exposure).
Q4  rule artifact Delta|b| at 0 ps is smaller in C0 than in the frozen stream in every tube, by a factor
    0.65-0.85 pooled (over-subtraction unchanged, tail signal x~4/3).
Q5  the pedestal law (slope and intercept frozen from Fig. 8c) predicts the eight C0 Delta|b| values
    from C0's own measured over-subtraction within 3 sigma each (sigma = bootstrap (+) 0.37 pp scatter).
Q6  retention Delta|b|(388, C2) / Delta|b|(0, C0), inverse-variance pooled over the eight tubes:
    0.42-0.54, i.e. higher than the frozen-baseline 0.359 and still below the registered 0.56-0.74.
Every prediction is scored and reported, pass or fail. Results within the frozen stream (mechanism,
CRT invariance, pedestal law, acceptance, rate replay) are paired on identical events and are not
re-derived; their magnitudes are reported with the yield caveat and the C0 values alongside.
