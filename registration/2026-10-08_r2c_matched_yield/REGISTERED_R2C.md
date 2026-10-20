# R2c: matched-yield, paired core-shell arms (review round 8, items 14 and 15)

Registered 2026-10-08, before any R2c data exist.

## Why
- Referee round 8, item 14: the stopping-point lookup selects both the lifetime and the o-Ps yield, so the
  production core-shell arms (R2b, 2026-08-13) do not isolate the lifetime; a matched-yield sensitivity is wanted.
- Item 15: a paired uncertainty of the arm difference is wanted. The R2b arms cannot be paired: each shard is a
  sequential Geant4 run, and the arms consume different numbers of random draws as soon as the source and
  stopping materials differ in I3 (first divergence at record 1657 of pair A and 71 of pair B in shard 1).

## Design
- Geometry, macros, seeds, shard layout, and primaries as R2b: r2_pairA.mac (core adipose_moskal2023 2.72 ns,
  shell brain_tumor_moskal2024 1.77 ns) and r2_pairB.mac (core adipose_avachat2024 2.54 ns, shell
  water_shibuya2020_n2 1.9215 ns), 10^7 primaries, 16 shards, seeds 1000000 + t.
- Materials: database v8 (pli_materials.json 23f607ea) with I3 = 0.2625 for every o-Ps-intensity entry, so that
  every material, including the total-formation structural entries (I3 = 0.35), has P_oPs = 0.2625,
  P_pPs = 0.0875, P_direct = 0.65. Lifetimes unchanged.
- PLI_DISABLE_3GAMMA=1, so that the decay-mode draw (P_3g = tau/142.05) cannot differ between arms.
- Arms: PLI_LIFETIME_MODE = transport and source for each pair; production binary e37be8b4.
- Readout as R2b: tail mean above 1 ns minus 1 ns on the true delays of the PLI_R2 records.

## Predictions
- P1 (pairing): within each pair the two arms produce the same number of records and identical
  (source material, annihilation material, channel) sequences, and every o-Ps delay ratio transport/source equals
  tau(annihilation material)/tau(source material) to the print precision.
- P2 (matched yield): Delta tau_hat = tau_hat(transport) - tau_hat(source) equals the first-order relation with
  equal yields, -82.06 ps (pair A) and -57.18 ps (pair B) (computed 2026-10-08 from the R2b transition matrices),
  within 3 ps each.
- P3 (paired precision): the paired standard error of Delta tau_hat is at least a factor 2 smaller than the
  independent combination of the two arms' standard errors.

## Not predicted
- The first-order relation evaluated on the R2c transition matrices is reported alongside, as a consistency check.
