# Registration records

The registration records of the analyses behind "Material-dependent positronium lifetimes in GATE and prompt-selection bias" (D. Dhar and J. Peter, submitted to Physics in Medicine and Biology). A registration record is a file in which a prediction, tolerance or falsifier was written down before the measurement it concerns. Every record here is a copy of the original with its modification time preserved. In 12 records, words of the project's internal working notes were replaced by neutral wording ("the authors", "decision", "request", "specified"). No prediction, criterion, number or date was changed. `NEUTRALIZED.tsv` lists these records with the sha256 of the original and of the copy, `DATES.tsv` keeps the sha256 of every original, and the originals are available from the authors on request.

## Contents

| path | what |
|---|---|
| `2026-08-20_R1_estimator/` ... `2026-10-08_r2c_matched_yield/` | the 25 records of the protocol (2026-08-20 onward): 23 Markdown files and 2 analysis scripts whose docstrings hold the registration |
| `pre_protocol/` | 6 records written 2026-08-15 to 2026-08-19, before the protocol started |
| `DATES.tsv` | modification time of each original, to the nanosecond, with size and sha256 |
| `NEUTRALIZED.tsv` | the 12 records whose wording was changed, with the sha256 of the original and of the copy |
| `SHA256SUMS` | sha256 of every other file in this folder (`sha256sum -c SHA256SUMS`) |

Folder names are those of the project results tree (`reconstruction/results/`), so each record keeps the path of the analysis it registers. The result files of those analyses are not in this folder. The subset the figures need is in `../figures/source_data/`, under the same folder names.

## 1. The protocol

**Start.** Pre-registration began on 2026-08-20 with `2026-08-20_R1_estimator/REGISTERED_R1.md`. From then on, before an analysis or simulation whose outcome was to be reported, a record was written into that analysis' result folder. It names the measurement, the expected value, band or direction, and what outcome would count as a failure. Outcomes were written into separate findings files, not into the record.

**Append-only.** A record is not edited once data it concerns exist. Five records received later sections, each headed with its own date and written before the data that section concerns: `REGISTERED_C.md` (per-configuration numbers "fixed 2026-08-24"), `REGISTERED_C2.md` (section "APPENDED 2026-08-25"), `REGISTERED_LADDER.md` (twelve-rung extension), `REGISTERED_ORACLE_GAP.md` (second falsifier "added 2026-08-24") and, before the protocol, `REGISTERED_N2_N3.md` (section N3b).

**Dates.** The project was not under version control. The filesystem modification time (mtime) of each file is its only date record. `DATES.tsv` lists the mtimes of the originals (local time, CEST, UTC+2). For a file with an appended section the mtime is the time of the last append, so its earlier sections are dated only by their own text. Copies that should keep the dates need `cp -p`, `rsync -t` or `tar` (git does not keep mtimes).

**Predictions 18 and 25** of section 2 are not separate files. They are the docstrings of the analysis scripts that test them, `2026-08-26_A2_cond2/sigma_geo_sim1.py` (prediction 18) and `2026-08-29_T6_diagnostics/t6_prior_test.py` (prediction 25, labelled P-T6-1, with P-T6-2 and P-T6-3 beside it). Each docstring says the prediction was registered before running, and in each case the script's output file is 6 s newer than the script. The scripts are shipped whole, unchanged. They import modules from the results tree by absolute path and are not meant to run from this folder.

**Paths.** `/data/dhritman/...` and `/home/dhritman/...` in the records are locations on the project workstation (bulk data and the project tree). They are not shipped and are not needed to read the records.

**Internal references in the records.**

- "decision", "request", "the authors'" (prediction, expectation): decisions, requests and predictions of the authors.
- "reviewer item N", "reviewer N", "review round N", "referee round N", "referee report": numbered points and rounds of internal pre-submission reviews of the manuscript, not the journal's peer review.
- "page N", "book page": pages of the project's internal results book, not shipped. "freeze": the project's internal append-only register of headline numbers, not shipped.
- Item labels such as V2, R1, G3, A1, B2, B5, C1, C2, C0, D1, D2 and N1 to N4 are the work items of the project log. "Sim-1" is the eight-tube phantom. "The frozen stream", "June" and "discovery stream" denote the eight-tube 0 ps acquisition of June 2026, which predates the correction of the o-Ps intensity convention of 2026-08-14. The paper calls it the earlier-yield acquisition (Supplementary Note S8). "C0" is the corrected-yield 0 ps rerun of 2026-09-30 and "C2" the 388 ps eight-tube arm.

**Pre-protocol records.** The six files in `pre_protocol/` were written in the same form between 2026-08-15 and 2026-08-19, before the protocol started. None is scored in the paper. `REGISTERED_V2.md` holds the water measurement from which the band of prediction 7 (0.56 to 0.74) was taken.

**Not shipped.** Three ledger documents of the results tree summarize the records at earlier dates (`2026-08-26_C_methods_facts/C2_holdout_statement.md`, `2026-08-27_B1_route2/RIDERS_R2_R3_R4.md`, `2026-08-27_item4/TABLE_predictions.md`). They are not registrations and are superseded by section 2 and by this index. The first version of the table among them differs from section 2 on predictions 3, 15, 21 and 23.

### Timing details

These facts come from the mtimes in `DATES.tsv` and of the result files in the same folders.

- **Prediction 9 (C2-b).** The two-factor retention model and C2-b were stated in section 2 of `REGISTERED_C2.md`, dated 2026-08-24 in its text, before the 388 ps eight-tube arm existed. The same section records that the second factor, the lever, did not reproduce the number it was to reproduce, and holds C2-b unscored. The C2 retentions were produced on 2026-08-26 at 08:04 (`2026-08-24_C2_sim1_388/c2_result.json`). The lever was then defined by a derivation of the first author dated 2026-08-26 (route 1, quoted in `2026-08-27_B1_route2/REGISTERED_ROUTE2.md`, the note itself is not in the results tree) and verified by the independent route registered in that file on 2026-08-27 at 00:55:29 (predictions 19 to 21). C2-b was scored four minutes later (`c2b_score.json`, 00:59:39, 10.1% RMS against the registered 15%). The model was therefore stated on 2026-08-24, before the 388 ps arm existed, and its second factor was defined on 2026-08-27, when the two routes agreed, after the retentions it is scored on had been produced. Section 2 scores prediction 9 as post hoc.
- `REGISTERED_LADDER.md` (mtime 2026-08-26 10:45:21) is newer than the six-rung result `ladder_result.txt` (01:28:56), because the twelve-rung section was appended at 10:45, before `mu_ladder.json` (11:01:06). Predictions P-L1 to P-L3 are dated by the text only.
- `REGISTERED_C2.md` (mtime 2026-08-26 00:19:19) carries the 2026-08-25 append. C2-a to C2-c are dated by the text of section 2 only.
- `REGISTERED_C.md` (mtime 2026-08-24 10:48:16, the 2026-08-24 append) precedes the C1 run script `2026-08-24_C1_nema388/run_c1.sh` (10:49:46).
- `REGISTERED_COND2.md` (12:26:17) is 104 s newer than `sigma_geo_sim1.json` (12:24:33), which it cites as the outcome of prediction 18, and older than the first P-C2 output (`cond2_component_shift.json`, 12:28:49).
- `REGISTERED_TANN.md` has the same mtime second as `tann_equivalence.py` (11:07:15). The output `tann_equivalence.json` followed at 11:07:49.
- `REGISTERED_INFLIGHT.md` (20:30:35) is newer than the physics-level output in its folder (`physics_level.json`, 20:29:07), which the record reports as done before it. Its predictions concern the GATE arm, whose comparison file `compare_inflight.json` is from 23:58. The record's header gives 20:35.
- `REGISTERED_C0.md` (21:52:43, header "21:50") is 70 s newer than the chain script `run_c0_chain.sh` (21:51:33). The comparison `c0_compare.json` is from 2026-10-01 03:28.
- `REGISTERED_C0_LADDER.md` states that it is not a registration in the sense of the protocol but a post hoc replication, and the paper counts it among the analyses that were not registered (Methods). Its expectations and outcomes are in note c of `docs/EIGHT_TUBE_LINEAGE.md`.

## 2. The scored predictions

The paper summarizes these predictions in Supplementary Note S2. Up to revision 2.20 of the manuscript they were Table S4 of the supplement, and the table below is that table with the record and label of each row. Rows 1 to 25 belong to the campaign of the rate replay of the earlier-yield stream, the realistic-timing runs and the estimators, rows 26 to 31 to the in-flight arm (Supplementary Note S8), rows 32 to 37 to the corrected-yield arm and rows 38 to 40 to the matched-yield core-shell rerun (Supplementary Note S5). pp: percentage points of the configured lifetime.

### Scoring rules

- A criterion is met only if every part of it is met, compared with the unrounded value, and split if some parts are met. A criterion stated for every tube, region or readout is one part, which is not met if any of them misses it.
- A correlation at the 5% level was considered significant for criteria of no correlation.
- Prediction 9 was stated before the data of the 388 ps arm existed, but its second factor was fixed only after the retention data it is scored on had been seen (section 1, timing details). It is therefore scored as post hoc.
- The registration of prediction 3 did not fix the tube set, which is scored on every tube with an estimate (0.989), whereas the six tubes common to all rungs, which give 1.009, became the tube set of prediction 4 and of the later replays.
- Predictions 15 to 17 were scored with the clipped pseudo-likelihood for which they were registered, on the bootstrap weights of the registered analysis (9.31 ± 0.40 pp and 1.95 times the truncated-mean shift). Supplementary Note S7 gives 9.20 ± 0.40 pp and 1.92 times for the same estimator on the 400 paired resamples it shares with the joint likelihood. The joint likelihood, at 10.46 ± 0.23 pp and 2.19 times, leaves every outcome unchanged.
- Predictions 12 and 14 are scored on the folded Δ\|b\| for which they were registered. On the signed difference prediction 12 would be met (largest paired-to-independent ratio 0.23) and prediction 14 would remain split (scatter 0.29 pp, 21% of the residual variance).
- The predictions were not corrected for multiplicity and serve as descriptive diagnostics, not as confirmatory tests.

Of the 24 prospective predictions of rows 1 to 25, 9 met their criterion, 10 did not, 4 split and 1 did not trigger. Of the 15 predictions of rows 26 to 40, 11 met their criterion, 3 did not and 1 split. Over all 40 rows, 20 met, 13 did not, 5 split, 1 did not trigger and 1 (prediction 9) is post hoc.

### Table of scored predictions

| # | prediction | registered | measured | outcome | record | label |
|---|---|---|---|---|---|---|
| 1 | rate-replay k = 1 rung reproduces the untouched extraction | byte-identical | byte-identical, both rules | met | `2026-08-26_rate_scan/REGISTERED_LADDER.md` | P-L1 |
| 2 | μ strictly proportional to k | proportional | μ/k constant, rel. sd 8×10^-10 | met | `2026-08-26_rate_scan/REGISTERED_LADDER.md` | P-L2 |
| 3 | replay exponent, six rungs | 1.0 to 1.5 (nearer linear than quadratic) | 0.989 (all tubes, five usable rungs), 1.009 (six common tubes), twelve rungs 0.999 (re-extraction 0.997, 95% interval 0.973 to 1.017) | not met (narrow) | `2026-08-26_rate_scan/REGISTERED_LADDER.md` | P-L3 |
| 4 | twelve-rung exponent against six-rung | within ±0.05 of 1.009 | 0.999 | met | `2026-08-26_rate_scan/REGISTERED_LADDER.md` (twelve-rung section) | P-L4 |
| 5 | k = 0.15 rung degenerate | excludable as bound | not degenerate, on the line | did not trigger | `2026-08-26_rate_scan/REGISTERED_LADDER.md` (twelve-rung section) | P-L5 |
| 6 | mixed-phantom exponent (μ² mechanism) | 2.0±0.5 | 1.42, confounded | not met | `2026-08-24_C2_sim1_388/REGISTERED_C2.md` (section appended 2026-08-25) | P-R3 |
| 7 | NEMA retention | 0.56 to 0.74 | 0.4655 (lifetime units) | not met | `2026-08-22_C_registration/REGISTERED_C.md` | P-C1 |
| 8 | retention monotone in τ0 | monotone | r = -0.34 | not met | `2026-08-24_C2_sim1_388/REGISTERED_C2.md` | C2-a |
| 9 | two-factor retention model | ≤15% RMS | 10.1% RMS | met (post hoc) | `2026-08-24_C2_sim1_388/REGISTERED_C2.md` | C2-b |
| 10 | one-factor model fails | ≥4 of 8 tubes | 8 of 8 outside | met | `2026-08-24_C2_sim1_388/REGISTERED_C2.md` | C2-c |
| 11 | inter-rule correlation | ρ > 0 everywhere | +0.81 to +0.99 (fifteen eight-tube and NEMA regions) | met | `2026-08-26_B2_paired/REGISTERED_B2.md` | P-B2-1 |
| 12 | paired/independent σ | < 1 everywhere | one region 1.29 (folded) | not met | `2026-08-26_B2_paired/REGISTERED_B2.md` | P-B2-2 |
| 13 | ρ falls when the direct component is displaced | falls | +0.964 against +0.960 | not met (null) | `2026-08-26_B2_paired/REGISTERED_B2.md` | P-B2-3 |
| 14 | systematic floor exists and dominates | both | 0.37 pp exists, 27%, not dominant | split | `2026-08-26_B2_paired/REGISTERED_B2.md` | P-B2-4 |
| 15 | component-estimator shift | 2 to 4 pp, positive | 9.31±0.40 pp | split (positive, magnitude not met) | `2026-08-26_A2_cond2/REGISTERED_COND2.md` | P-C2-1 |
| 16 | shift resolved | \|z\|>3 | z = 23 | met | `2026-08-26_A2_cond2/REGISTERED_COND2.md` | P-C2-2 |
| 17 | component shift against truncated mean | smaller | 1.95× larger | not met | `2026-08-26_A2_cond2/REGISTERED_COND2.md` | P-C2-3 |
| 18 | eight-tube σ_geo against NEMA | smaller | 353.2 against 351.4 ps | not met | `2026-08-26_A2_cond2/sigma_geo_sim1.py` (docstring) | none (sigma_geo of Sim-1 against NEMA) |
| 19 | pedestal-response factor of prediction 9, second implementation, pure points | agree to 1% | 0.004%, 0.028% | met | `2026-08-27_B1_route2/REGISTERED_ROUTE2.md` | P-R2-1 |
| 20 | pedestal-response factor of prediction 9, second implementation, mixture | agree to 1% | 0.000% | met | `2026-08-27_B1_route2/REGISTERED_ROUTE2.md` | P-R2-2 |
| 21 | 0.1996 recoverable as a wrong step | recoverable as a specific wrong step. Leading hypothesis: configured instead of observed lifetime, or uninverted window mean | recovered, neither hypothesis reproduces it | split | `2026-08-27_B1_route2/REGISTERED_ROUTE2.md` | P-R2-3 |
| 22 | simplified t_ann equivalent | \|δ\|<2σ all VOIs | 3 of 8 beyond, pooled z = -4.7 | not met | `2026-08-27_tann/REGISTERED_TANN.md` | P-T1 |
| 23 | per-triple difference of the two delay variables: spread and mean | spread of order σ_geo (300 to 360 ps), mean +10 to +40 ps | spread 359.2 ps (σ_geo 353.2 ps), mean -35.1 ps | split (spread met, mean of opposite sign) | `2026-08-27_tann/REGISTERED_TANN.md` | P-T2 |
| 24 | δ against τ_cfg | no correlation | r = 0.39, p = 0.33 (8 regions) | met | `2026-08-27_tann/REGISTERED_TANN.md` | P-T3 |
| 25 | removing the intensity prior halves the healthy-brain deviation | from 6.26% to at most 3.2% in absolute value | -6.30%, no tube moves > 0.08 pp | not met | `2026-08-29_T6_diagnostics/t6_prior_test.py` (docstring) | P-T6-1 |
| 26 | 511-window singles | -3 to -1% | -1.81% | met | `2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md` | P1 |
| 27 | prompt-window singles | within ±0.5% | +0.14% | met | `2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md` | P2 |
| 28 | event-identity signal triples per tube | -3.5 to -1.5% | mean -2.33%, one tube -1.46% | not met (narrow) | `2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md` | P3 |
| 29 | lifetimes (event identity, all-prompts, closest-prompt) | \|z\|<2 each | z = -1.28, +1.98, +2.07 | not met (narrow) | `2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md` | P4 |
| 30 | Δ\|b\|, pooled | \|z\|<2 | z = +0.4 | met | `2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md` | P5 |
| 31 | accidental-prompt rate μ | within ±0.5% | +0.14% | met | `2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md` | P6 |
| 32 | signal ratio in every tube | 1.25 to 1.42 | 1.29 to 1.32 | met | `2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md` | Q1 |
| 33 | event-identity lifetimes unchanged | \|z\|<3 per tube, < 2 pooled | pooled z = -5.8, largest \|z\|=4.7 (within 75 mm of the decay: pooled +0.3, largest 1.3) | not met | `2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md` | Q2 |
| 34 | μ against the realistic-timing arm | within ±0.5% | -0.01% | met | `2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md` | Q3 |
| 35 | Δ\|b\| smaller in every tube, pooled factor | 0.65 to 0.85 | seven of eight, pooled 0.79 | split | `2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md` | Q4 |
| 36 | pedestal law in every tube (0.37 pp scatter) | within 3σ | largest pull 2.0 | met | `2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md` | Q5 |
| 37 | pooled retention | 0.42 to 0.54 | 0.437 | met | `2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md` | Q6 |
| 38 | arms matched event for event | identical sequences | identical in both pairs | met | `2026-10-08_r2c_matched_yield/REGISTERED_R2C.md` | P1 |
| 39 | Δτ̂ of pairs A and B | within 3 ps of -82.1 and -57.2 ps | -82.2 and -56.8 ps | met (A and B) | `2026-10-08_r2c_matched_yield/REGISTERED_R2C.md` | P2 |
| 40 | paired against independent σ | at least two times smaller | factors 6 and 9 | met | `2026-10-08_r2c_matched_yield/REGISTERED_R2C.md` | P3 |

The record labels P-C2 (in `REGISTERED_C.md`) and P-C2-1 to P-C2-3 (in `REGISTERED_COND2.md`) are unrelated, as are P-R2 (in `REGISTERED_C2.md`) and P-R2-1 to P-R2-3 (in `REGISTERED_ROUTE2.md`).

### The corrected-yield replay

| paper | record | labels |
|---|---|---|
| not registered (Methods), Supplementary Note S7 and note c of `docs/EIGHT_TUBE_LINEAGE.md` | `2026-10-08_c0_ladder/REGISTERED_C0_LADDER.md` | E1 to E4 |

## 3. Index of every prediction

Every prediction, tolerance and falsifier in every record, with its outcome. For the predictions of section 2, the outcome is that of section 2, or for E1 to E4 that of note c of `docs/EIGHT_TUBE_LINEAGE.md`. For the others, the outcome was taken from the findings documents and result files of the same result folder and from the project log (not shipped), and is scored here with the rules of section 2: unrounded values against the registered criterion, a criterion with several parts is met only if every part is, and split if some are, and a criterion stated for every tube, region or readout is one part. Where a record gives a direction without a number, the measured value is quoted. Where a record or findings document scores a prediction itself, that score is quoted, and any difference from a literal reading of the criterion is stated.

Outcome categories:

- met, not met, split: as in the paper
- did not trigger: a conditional prediction or falsifier whose condition did not occur
- partly scored: some parts were measured and scored, others never measured
- not scored: no outcome was ever computed, with the reason where known
- no directional prediction: the record expects a quantity "may change" without a direction, the measured value is given
- withdrawn: withdrawn in the record itself before the measurement

### Summary

| | predictions | met | not met | split | did not trigger | partly scored | not scored | no direction | withdrawn |
|---|---|---|---|---|---|---|---|---|---|
| scored in section 2, rows 1 to 25 | 25 | 10 | 10 | 4 | 1 | | | | |
| scored in section 2, rows 26 to 40 | 15 | 11 | 3 | 1 | | | | | |
| corrected-yield replay (`docs/EIGHT_TUBE_LINEAGE.md`, note c) | 4 | 4 | | | | | | | |
| not in the paper, protocol records | 48 | 21 | 13 | 1 | 3 | 2 | 8 | | |
| not in the paper, pre-protocol records | 23 | 14 | 1 | 1 | | | 4 | 2 | 1 |
| total | 115 | 60 | 27 | 7 | 4 | 2 | 12 | 2 | 1 |

The met count of 10 for rows 1 to 25 includes prediction 9, which section 2 scores as post hoc. The 48 protocol predictions not in the paper are in 15 records: 12 records with no row in the paper (R1, G3, density, G6, residual structure, A1b, A6, B2 and B2b of 2026-08-22, ACC, ORACLE_GAP, P5) and three records whose other predictions are in section 2 (P-C2 and P-C3 in `REGISTERED_C.md`, P-R1 and P-R2 in `REGISTERED_C2.md`, P-T6-2 and P-T6-3 in `t6_prior_test.py`).

### Protocol records, 2026-08-20 onward

**`2026-08-20_R1_estimator/REGISTERED_R1.md`** (estimator swap on identical events). Findings: `findings_R1.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| R1-1 | the windowed, truncation-corrected tube estimator and the raw oracle estimator will not agree, and the question is whether they differ by less than the 9-16 ps offset | partly scored | 388 ps arm: they differ by -244 to -250 ps in three regions. Production arm: -15.1 ± 42.9, +5.9 ± 18.2, +28.4 ± 10.9 ps, which findings_R1.md calls underpowered for a 9-16 ps question |
| R1-2 | if the tube estimator reproduces the raw oracle within a few ps, the tube closure and the reconstruction offset are separate findings | not scored | condition not established (underpowered on the production arm) |
| R1-3 | if the tube estimator differs by about the offset, the D2 trend (-19.11 ± 6.27 ps/ns) should follow from propagating a constant per-event delay through the inversion | not met | a constant delay (and a symmetric smear) propagated through the inversion gives +0.22 ps/ns. The estimator absorbs additive delays exactly (findings_R1.md) |

**`2026-08-22_G3_consequences/REGISTERED_G3.md`** (consequences of paired errors). Results: `g3a_nema.json`, `g3b_pedestal.json`, `g3c_paired_estimator.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| G3-a | NEMA paired/independent sigma ratio near 0.14-0.19, significance 24-33 sigma | not met (narrow) | ratio 0.135, z = 34.3 (g3a_nema.json) |
| G3-b | pedestal-law slope stays near 0.99, chi2/dof rises well above 1 | met | slope 1.015, chi2/dof 7.31 over 20 regions, 3 beyond 3 sigma (g3b_pedestal.json). After the 2026-08-27 correction of the error statistic the project log gives chi2/dof 2.71 |
| G3-c | pairing the two estimators shrinks sigma substantially, and if the paired sigma falls below about 4 ps the D1 offset becomes testable | not met | paired/independent 0.72 to 0.81, paired sigma 10.9 to 40.6 ps (g3c_paired_estimator.json). The next record (`REGISTERED_G6.md`) describes this as shrinking "only to 0.72-0.81". The 4 ps condition did not occur |

**`2026-08-22_density_fix/REGISTERED_density.md`** (XAD-4 densities corrected). The folder holds the singles counts of both densities. The lifetimes are in the project's results book, page 20.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| D-1 | the extracted lifetime is unchanged within counting error | met | z = -1.11, -1.05, -0.76, +0.93 for the four samples (results book page 20) |
| D-2 | singles yield falls, the scattered fraction from within the tube falls and the unattenuated 511 keV flux rises (direction only) | partly scored | kept singles -1.94% (22,628 to 22,189, `singles_oldrho.json`, `singles.json`), 511 keV window -1.69%. The scattered fraction and the unattenuated flux were not measured separately |
| D-3 | if the lifetime moves by more than counting error, density is coupled to the lifetime | did not trigger | see D-1 |

**`2026-08-22_G6_channel/REGISTERED_G6.md`** (channel composition above the censor). Result: `g6_channels.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| G6-1 | the fit-window population is overwhelmingly o-Ps, p-Ps at the 1e-4 level, so no channel contaminant explains the offset | met | p-Ps fraction 1.93e-4, direct 0, shift -0.33 ps on 5,563,082 annihilations above 1 ns |
| G6-2 | composition is the same on and off axis within counting error | not scored | `g6_channels.json` holds one population, no on/off-axis split |

**`2026-08-22_residual_structure/REGISTERED.md`** (residual structure of the pedestal law). Findings: `findings.md`, result `residual_structure.json`. The findings' verdict: all three hypotheses fail to explain the residual structure.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| H1 | pulls sort by lifetime, hepatic and polycarbonate at opposite extremes | not met | the two extremes are 0.030 ns apart, r(tau, pull) = -0.060 (falsified in the record itself before the test) |
| H2 | abs(residual) rises monotonically with the randoms-fraction difference, outliers at the largest values | not met | r = +0.031, outlier ranks 20, 2 and 4 of 20 |
| P3a | residual correlates with (m_ped* - 2.0) with the sign the law dictates | met | r = +0.657 over 20 regions, predicted sign. The findings attribute it to noise shared by the two quantities (r = +0.280 over the 14 well-measured regions) |
| P3b | substituting the measured pedestal moment lowers the RMSE about the identity line below 3.74% | met | RMSE 3.74% to 2.66% |
| P3c | the 3-sigma outliers are the regions with the largest abs(m_ped* - 2.0) | not met | their ranks are 5, 9 and 6 of 20 (T8 and T4 have equal deviations and residuals of -0.91% and -7.88%) |

**`2026-08-22_residual_structure/REGISTERED_A1b.md`** (edge deficit of the delayed window). Findings: `findings.md` (addendum), `a1_edge_p3.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| A1b-P1 | the deficit appears at both window edges | met | 29% at the low edge, 23% at the high edge |
| A1b-P2 | falsifier: a deficit at the high edge only refutes the mechanism | did not trigger | see A1b-P1 |
| A1b-P3 | the edge roll-off does not correlate with the region's residual | met | r = +0.281 over 14 regions (p about 0.33) |

**`2026-08-22_A6_moskal/REGISTERED_A6.md`** (oxygen response against Moskal arXiv:2605.00021). Result: `a6_moskal_compare.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| A6-P1 | water: our lifetime change lower than his 5 ps, of order 3-4 ps | met | 3.17 ps |
| A6-P2 | adipose: ours several times smaller than his 48 ps | met | 6.38 ps, a factor 7.5 |
| A6-P3 | the ratio of the two changes tracks the ratio of the quenching increments within the tau^2 factor | met | 1.58 against 1.43 (water), 7.52 against 7.14 (adipose), within 11% |

**`2026-08-22_B2_estimator/REGISTERED_B2.md`** (synthetic gate of the three-component estimator). Findings: `findings_B2.md`, `b2_gate.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| B2-T1 | lifetime within 1% of truth at every grid point, N = 1e6 | met | worst 0.74% |
| B2-T2 | within 3 sigma at every grid point, N = 1e6 | met | worst 2.11 sigma |
| B2-T3 | o-Ps intensity within 0.02 of truth | met | worst 0.002 |
| B2-T2b | within 3 sigma at N = 1e5 | met | worst 1.34 sigma |

**`2026-08-22_B2_estimator/REGISTERED_B2b.md`** (the gated estimator on simulated water spectra). Findings: `findings_B2.md`, `b2b_real.json`, `b2b_range_scan.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| B2b-P1 | on the 388 ps arm the fit recovers 1.75-1.93 ns under both rules | not met | 1.347 ns (closest-prompt), 1.391 ns (all-prompts). The result moves by 407-625 ps with the end of the fit range |
| B2b-P2 | falsifier: a fit near 1.15 ns refutes the claim that the deficit belongs to the simple estimator | did not trigger as worded | the fit returned 1.35-1.39 ns. The findings nonetheless withdraw that claim because P1 failed |
| B2b-P3 | on the 0 ps arm the fit agrees with the truncated mean within its uncertainty | met | agreement to 37 and 19 ps (findings_B2.md scores it as held) |

**`2026-08-22_C_registration/REGISTERED_C.md`** (the 388 ps chain, band from V2 water). Results: `2026-08-24_C1_nema388/c1_result.json`, `2026-08-24_C2_sim1_388/c2_result.json` and `out_*/per_voi.csv`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| P-C1 | NEMA retains 0.56-0.74 of its 0 ps artefact, shift +0.039 to +0.052 ns | not met (row 7) | 0.4655 |
| P-C2 | each Sim-1 tube retains 0.56-0.74 of its own 0 ps DiD (per-tube bands in the record) | not met | retention 0.303 to 0.487, 0 of 8 tubes inside (c2_result.json). Its complement is row 10 |
| P-C3 | the truncated-mean lifetime falls by 0.4-0.7 ns at 388 ps in NEMA and in Sim-1 | split | NEMA: 528.8 ps (closest-prompt), 566.5 ps (all-prompts), met. Sim-1: 0.405 to 0.940 ns (closest-prompt) and 0.505 to 1.196 ns (all-prompts) against the June 0 ps stream, 10 of 16 tube-rule values above 0.7 ns, not met. The project log scores only the NEMA half. The Sim-1 values were computed for this index from `out_*/per_voi.csv` and the June biases in `2026-09-30_feedback3/c0_sim1/c0_compare.json` |

**`2026-08-23_acceptance/REGISTERED_ACC.md`** (acceptance term, re-fit, D2 test). Findings: `findings_ACC.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| ACC-R1 | the analytic acceptance (sigma_e = 64.3 ps plus a small geometric term) agrees with the measured edge shape | not met as registered | measured 290.6 ps, the 64.3 ps form misses the edge bin by 51 sigma. findings_ACC.md scores a repaired route that treats the geometric term as dominant (312.4 ps a priori, agreement 7%) as a pass |
| ACC-gate | the four B2 gate conditions hold on spectra generated with the acceptance | met | worst 0.69%, 1.04 sigma, o-Ps share within 0.001, 1.43 sigma at 1e5 |
| ACC-R2 | the re-fit returns 1.75-1.93 ns on the 388 ps water arm under both rules | not met | closest-prompt 1.720-1.754 ns (1 of 6 fit endpoints inside), all-prompts 1.747-1.816 ns (5 of 6). The findings call it partial |
| ACC-R3 | the endpoint swing falls below 3 times the statistical error | not met | closest-prompt 34.4 ps against 37.5 ps (met), all-prompts 68.8 ps against 39.5 ps (not met) |
| ACC-R4 | nothing is retuned per sample | met | one acceptance, one model for every arm |
| ACC-R5 | the propagated acceptance reproduces the D2 slope -19.11 ± 6.27 ps/ns (falsifier R7: off by more than the error) | not met | -136.0 ps/ns, 18.7 sigma. findings_ACC.md judges the test inapplicable because the D2 tube chain has no ±5 ns prompt box |
| ACC-R6 | the zero crossing near 2.01 ns is a consistency check on the +9 to +16 ps offset | not scored | not reported in findings_ACC.md |
| ACC-R8 | the edge deficit does not explain the pull outliers | met | T8 has the larger deficit and one ninth of T4's residual |

**`2026-08-23_B3/REGISTERED_ORACLE_GAP.md`** (the oracle sits 280-341 ps above the mixture prediction on the Mercolli tubes). The record registers an observation and two falsifiers for any mechanism proposed later.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| OG-F1 | a proposed mechanism must reproduce the near-constant gap across the four tubes | not scored | no mechanism was proposed and tested afterward |
| OG-F2 | it must also explain why the same algebra closed to within 0.4% on V2 water | not scored | as OG-F1 |

**`2026-08-24_P5_overflow/REGISTERED_P5.md`** (overflow-bin prompt selection). Findings: `findings_P5.md`, which retracts the measurement the same day.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| P5-a | the genuine-triple fraction falls by 5-15 percentage points | not scored | measurement retracted: the stored singles contain no photons between 562 and 1080 keV, so the overflow window could not be applied, and the original ROOT files were deleted |
| P5-b | the truncated-mean lifetime falls to 1.03-1.12 ns | not scored | as P5-a |
| P5-c | sigma within a factor 0.9-1.3 of the C1 value | not scored | as P5-a |

**`2026-08-24_C2_sim1_388/REGISTERED_C2.md`** (eight tubes at 388 ps).

| label | registered | outcome | measured (source) |
|---|---|---|---|
| C2-a | retention falls monotonically with tau_0 | not met (row 8) | r = -0.34 |
| C2-b | two-factor model within 15% RMS | met, post hoc (row 9) | 10.1% RMS |
| C2-c | one-factor model fails on at least 4 of 8 tubes | met (row 10) | 8 of 8 outside |
| P-R1 | the corrected run reproduces mu = 0.029 within ±15% | met | mu = 0.02924 on the delivered data (project log) |
| P-R2 | per-tube DiD returns inside the bands registered on 24 August | not met | DiD 3.44 to 5.98 pp, retention 0.303 to 0.487, 0 of 8 inside (c2_result.json) |
| P-R3 | the artefact scales as mu^n with n = 2 ± 0.5 | not met (row 6) | 1.42, confounded |

**`2026-08-26_rate_scan/REGISTERED_LADDER.md`**: P-L1 met, P-L2 met, P-L3 not met (narrow), P-L4 met, P-L5 did not trigger (rows 1 to 5).

**`2026-08-26_B2_paired/REGISTERED_B2.md`**: P-B2-1 met, P-B2-2 not met, P-B2-3 not met (null), P-B2-4 split (rows 11 to 14).

**`2026-08-26_A2_cond2/sigma_geo_sim1.py`** (docstring): sigma_geo of Sim-1 smaller than NEMA's 351.4 ps, not met, 353.2 ps (row 18).

**`2026-08-26_A2_cond2/REGISTERED_COND2.md`**: P-C2-1 split, P-C2-2 met, P-C2-3 not met (rows 15 to 17).

**`2026-08-27_B1_route2/REGISTERED_ROUTE2.md`**: P-R2-1 met, P-R2-2 met, P-R2-3 split (rows 19 to 21).

**`2026-08-27_tann/REGISTERED_TANN.md`**: P-T1 not met, P-T2 split, P-T3 met (rows 22 to 24).

**`2026-08-29_T6_diagnostics/t6_prior_test.py`** (docstring). Findings: `FINDINGS_T6.md`, `t6_prior_test.json`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| P-T6-1 | removing the Dirichlet prior shrinks abs(dev) on T6 by more than half | not met (row 25) | -6.26% to -6.30% |
| P-T6-2 | the shrink is largest on T6 of the eight tubes | met | abs(dev) grows on every tube without the prior. T6's change (-0.6%) is the smallest, which meets the criterion literally. FINDINGS_T6.md calls this pass vacuous once P-T6-1 failed |
| P-T6-3 | the fitted o-Ps intensity on T6 rises above its prior-constrained value | met | 0.296698 to 0.296737. FINDINGS_T6.md calls this pass vacuous |

**`2026-09-30_feedback3/inflight/REGISTERED_INFLIGHT.md`** (Supplementary Note S8): P1 met, P2 met, P3 not met (narrow), P4 not met (narrow), P5 met, P6 met.

**`2026-09-30_feedback3/c0_sim1/REGISTERED_C0.md`** (corrected-yield arm): Q1 met, Q2 not met, Q3 met, Q4 split, Q5 met, Q6 met.

**`2026-10-08_c0_ladder/REGISTERED_C0_LADDER.md`** (note c of `docs/EIGHT_TUBE_LINEAGE.md`, self-declared post hoc replication): E1, E2, E3, E4 met.

**`2026-10-08_r2c_matched_yield/REGISTERED_R2C.md`** (Supplementary Note S5): P1 met, P2 met (both pairs), P3 met.

### Pre-protocol records, 2026-08-15 to 2026-08-19

**`pre_protocol/2026-08-13_overnight/R6_brain_tumor_rerun/REGISTERED_PREDICTION.md`** (brain-tumor re-run).

| label | registered | outcome | measured (source) |
|---|---|---|---|
| R6-0 | +33% o-Ps yield lowers r_min | withdrawn | withdrawn in the record before the run (the two brain materials were not affected by the convention fix) |
| R6-1 | detection sigma and r_min(3 sigma) reproduce the June values within about 1% (9.693 mm at 302.7 ps, 8.913 mm at 228.1 ps) | not scored | the folder holds the per-radius detection sigma (`recovery_points.csv`: 2.12 sigma at 8 mm and 5.00 at 12 mm for 302.7 ps, 2.94 and 5.76 for 228.1 ps) and a recovery-curve r_min at recovery fraction 0.5, a different quantity. No r_min(3 sigma) and no comparison with the June values is on record |

**`pre_protocol/2026-08-13_overnight/V1_vipman_postcorrection/REGISTERED_PREDICTION.md`** (VIP-Man after the convention fix). Outcomes from the project log (entries of 2026-08-15/16), `B_voxel_sweep.csv` and `2026-08-13_overnight/V2_migration/findings.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| V1-1 | Table 3 residuals reproduce within the seed spread (0.80% mean, 1.98% worst) | met | mean abs change 0.630% (project log). The worst tissue is not stated |
| V1-2 | the migration budget is unaffected (22.6% within 1 mm, per-material net fluxes) | met | 22.58 ± 0.013%, all 8 net fluxes within 0.7 pp |
| V1-3a | activity reconstruction unchanged | not scored | no pre/post comparison on record |
| V1-3b | per-voxel lifetime maps unchanged | met | as V1-1 |
| V1-3c | blind region lifetimes 1.80-2.45 ns unchanged | not scored | not re-derived after the fix |
| V1-3d | K = 2 model-order selection may change | no directional prediction | K unchanged (K = 2 passes 32 of 32 seeds before and after) |
| V1-3e | component lifetimes E0, E1 may change | no directional prediction | E0 1.5057 to 1.4943 ns, E1 2.8768 to 2.9130 ns |
| V1-3f | AUCs survive | not scored | blocked: the per-batch oracle log needed to score them was deleted by the run recipe |

**`pre_protocol/2026-08-17_derivation_test/REGISTERED_N2_N3.md`** (over-subtraction). Findings: `findings_N2_N3_N4.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| N2-1 | all-prompts: subtracted and true accidentals agree within counting error | met | per-ROI differences -4.7 to +2.1 sd (the -4.7 sd in the disk phantom), which the findings call consistent with zero |
| N2-2 | closest-prompt: the subtraction exceeds the true accidentals | met | whole phantoms (background regions): +18.1 to +68.4 sd. Per region: +8.3 to +68.4 sd in every region of the eight tubes (discovery stream), the 20 mm water cylinder and the disks, and -1.5 to +3.6 sd in the NEMA spheres and central background, where one region, the 13 mm sphere (62 triples), lies below zero (-7 counts, -1.5 sd). The record states the prediction per phantom, read per region it would be split (`oversub_*_closest.json`) |
| N3 | abs(bias) is one monotonic function of phi, all four phantoms on one curve | split | one curve on the closest-prompt arm (R^2 0.926, between/within F = 0.03, NEMA out of sample at -0.06 sigma). Out of sample on the all-prompts arm it under-predicts by a median factor 2.5, and within Sim-1 phi is collinear with I3 (rho = -0.898), a failure mode the record names |
| N3b | the rule difference of abs(bias) is one monotonic function of the rule difference of phi, common to all phantoms, through the origin | met | slope 114.7 through the origin, R^2 0.980, per-phantom slopes 113.8, 127.4, 125.2 (free intercept +0.67%) |

**`pre_protocol/2026-08-19_V2_earlytime/REGISTERED_V2.md`** (direct annihilation at 388 ps). Findings: `findings_V2.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| V2-1 | absolute bias worsens under both rules | met | bg region: closest-prompt -0.386 to -0.696 ns, all-prompts -0.297 to -0.647 ns, against each arm's own oracle |
| V2-2 | the rule artefact shrinks and does not vanish | met | in all four regions. The stated basis (growth of the genuine count) predicts 71-76% retention against 56-74% measured, which the findings call quantitatively incomplete |
| V2-3 | the artefact retains more than half its magnitude | met | 56.1% to 73.7% |

**`pre_protocol/2026-08-19_V3_oxygen/REGISTERED_V3.md`** (oxygen sweep). Findings: `findings_V3.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| V3-1 | every point lies on the analytic line within its own error | not met | common offset +15.62 ± 1.66 ps (9.4 sigma), no slope |
| V3-2 | the seven physiological points do not separate | met | as registered |
| V3-3 | the three extended points separate cleanly | met | z = -2.24, -11.4, -33.1, 207 ps measured against 219 ps predicted at 3000 mmHg |

**`pre_protocol/2026-08-19_V3b_offset/REGISTERED_V3b.md`** (where the offset enters). Findings: `findings_V3b.md`.

| label | registered | outcome | measured (source) |
|---|---|---|---|
| V3b-1 | the sampled lifetime equals 1.840 ns within its error | met | -0.10 ± 0.78 ps |
| V3b-2 | the reconstructed oracle returns about 1.856 ns, reproducing the offset | met | +11.94 ± 5.46 ps in the largest region (+26.6 ± 21.4 and -1.0 ± 9.2 ps in the two smaller ones, no tolerance registered) |
| V3b-3 | the offset is therefore added by the reconstruction | met | the sampler is exact in every region and the offset enters downstream |
