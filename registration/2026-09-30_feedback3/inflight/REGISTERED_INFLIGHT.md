# In-flight annihilation sensitivity arm, registered BEFORE the GATE run (2026-09-30 20:35)

Feedback round 3 asks whether suppressing Geant4 in-flight annihilation changes any headline
observable. Hybrid physics: shipped PLI list, but the stock G4eplusAnnihilation is KEPT for its
PostStep (in-flight) DoIt and only its AtRest DoIt is deactivated (PLI_INFLIGHT_STOCK=1, separate
binary software/install/gate-inflight, sha256 cc552a68...). Gate: with the switch OFF the new binary
reproduces the production binary (e37be8b4...) on a 0.004 s Sim-1 probe, every numeric Singles and
Coincidences branch identical (gate_probe/probe_compare.json).

Physics level (pli_inflight_test, 1e6 positrons, done before this registration):
44Sc in-flight 2.265 +- 0.015 % hybrid vs 2.244 +- 0.015 % stock (mean KE 420 keV);
124I 2.896 % vs 2.887 % (520 keV). Delayed-tail (>1 ns) yield hybrid/shipped 0.972 +- 0.003 (44Sc).
Only 2.5 % of in-flight pairs have both photons in [460,562] keV.

Arm: identical to C2 (2026-08-24_C2_sim1_388): same macro, 388 ps direct delay, time-driven,
10 batches x 12 shards + 1 x 8, 0.137933 s per shard, same seeds, same converter windows.
Baseline = C2's own singles (/data/dhritman/c2_sim1_388). Only intended difference: in-flight.

Predictions (hybrid minus C2):
P1  in_511-flagged singles fall by 1 to 3 % relative (in-flight photons leave the 511 window).
P2  prompt-window singles change by less than 0.5 % relative.
P3  genuine (truth-tagged) triples per tube fall by 1.5 to 3.5 % relative.
P4  per-tube lifetime (truth-tagged, all-prompts, closest-prompt readouts): pooled difference
    consistent with zero (|z| < 2); per-tube differences consistent with the independent-arm noise.
P5  rule artifact Delta|b| = |b_closest| - |b_all|: pooled difference consistent with zero (|z| < 2).
P6  prompt-competition rate mu changes by less than 0.5 % relative.
Every prediction is scored and reported, pass or fail.
