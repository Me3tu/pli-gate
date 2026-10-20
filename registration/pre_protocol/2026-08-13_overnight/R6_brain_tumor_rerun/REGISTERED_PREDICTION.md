# Registered prediction — brain-tumor re-run

**Written 2026-08-15, BEFORE the run is launched.** The run starts only after
VIP-Man completes (~19:45). Nothing below may be edited after launch; the
measured outcome goes in a separate section of the findings.

---

## Withdrawn prediction

I previously wrote that the re-run should show

> "+33 % o-Ps yield → σ_f falls by √1.333 → detection σ rises ~15.5 % → r_min falls"

**This is wrong and is withdrawn.** It applied the R2 convention-fix yield gain
to a phantom that the fix does not touch.

`gate/geometry/phantom_brain_tumor_visibility.mac` sets exactly two materials:

* `brain_healthy_moskal2024` — the phantom body
* `brain_tumor_moskal2024` — all five inserts (r = 4, 8, 12, 16, 20 mm)

Both were **already declared `ops_intensity`** before the correction, and both
appear in the *untouched* group of the Sim-1 pre/post table, where they were
**measured**, not assumed:

| material | o-Ps yield pre → post | Δ |
|---|---|---|
| `brain_healthy_moskal2024` | 20.59 % → 20.60 % | **+0.0 %** |
| `brain_tumor_moskal2024` | 17.08 % → 17.07 % | **−0.1 %** |

The convention fix moved three tubes by ~+32 % and left five alone. Both brain
materials are in the second group. I had the refuting data in a table I wrote
myself and did not apply it.

## The actual prediction

The only mechanism that should move these numbers is the **three-photon
channel**, which is new in the frozen binary. From the R1 measurements:

| input | value | source |
|---|---|---|
| P(3γ \| o-Ps) | 1.5741 % | R1 production measurement |
| 3γ events retaining two photons in [460, 562] keV | 5.70 % | R1 |
| **o-Ps events removed from the 2γ sample** | **1.4844 %** | 0.015741 × (1 − 0.057) |

Propagating through the detection model `σ_f = C/√N_ROI`, `detection σ = f/σ_f`,
`σ_f ∝ r^(−3/2)`:

| quantity | predicted change |
|---|---|
| N_ROI | **−1.48 %** |
| σ_f | **+0.75 %** |
| detection σ | **−0.75 %** |
| r_min(3σ) | **+0.50 %** |

Absolute r_min targets:

| CRT | June value | predicted |
|---|---|---|
| 302.7 ps | 9.645 mm | **9.693 mm** |
| 228.1 ps | 8.869 mm | **8.913 mm** |

## Pass criterion

**Detection σ and r_min reproduce the June values to within ~1 %.** A result
that reproduces June is the *expected and correct* outcome and will be reported
as confirmation, not as a null or a failure.

**A movement larger than ~1 % is the finding**, and would mean one of:

* a material in the phantom is convention-affected after all (contradicting the
  measured Sim-1 rows above);
* the 3γ channel interacts with the ROI/shell estimator more strongly than a
  flat 1.48 % count loss;
* an exposure or configuration difference between the runs that has not been
  identified.

In that case: stop, report, and do not fold the number into the manuscript
until the cause is established.

## Why the run happens at all

**Provenance, not physics.** The purpose of the exercise is that one frozen
binary produces every number in the paper. Accepting the analytic projection
here would reintroduce the mixed provenance that took two days to remove, and
would leave the only inferred number in an otherwise measured results section.
A 5 h run that lands ~01:00 and blocks nothing is a cheap price for that.

## Configuration

Identical to the R5 brain-tumor run except that `singles=True` is now set, so
the ROOT is converted before deletion and `M4_recovery_curve.py --singles` can
run.

* macro: `2026-08-13_overnight/macros/brain_tumor_te_poisson_1B_truetime.mac`
  (archived-path fix; sub-macro byte-identical, sha256 `c29860a704bad160…`)
* **16 batches × 61.2 s** — matching `run_summary.json` batch_ids 0–15, *not*
  the 17 I ran in R5 from a mis-counted directory listing
* binary `3522d567…`, base database `2a26b68e…`
