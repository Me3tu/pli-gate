# VIP-Man post-correction — registered predictions

Written **before** the voxel sweep, migration check and blind re-derivation run.

---

## 0. Two points from the per-tissue result, recorded explicitly

### 0.1 The purity factors disagree in the RIGHT direction

| phantom | implied tail purity p |
|---|---|
| Sim-1, five already-declared tubes | **0.938** |
| VIP-Man, eight tissues | **0.923** |

This is not merely "agreement to 1.5 pp". **VIP-Man should have the lower tail
purity**, and does. The anthropomorphic phantom carries scatter and
cross-boundary positron migration that tubes suspended in air do not, so a
larger fraction of its tau > 1 ns window is I3-independent background. A purity
that came out *higher* for VIP-Man than for the tube phantom would have been the
result to worry about — it would have implied the dilution model was fitting
something other than contamination.

### 0.2 The scale check is load-bearing, not a footnote

`simple_recon` computes `m = w / max(x, 1e-3)` — lifetime-weighted image over
activity image — i.e. the mean delay over **all** annihilations, not over o-Ps
events only. Hence

    m  ~=  P(o-Ps) x tau_oPs  ~=  0.20 x 2 ns  =  0.4 ns

and the observed per-tissue m is 0.30-0.39 ns (May) rising to 0.40-0.51 ns (O4).
The scale matches, which is what **identifies the quantity**.

Without that identification, a +31 % shift in a per-tissue reconstructed
quantity is alarming — it reads as a lifetime catastrophe. With it, m is
proportional to o-Ps yield at fixed lifetime and the shift is exactly the
expected convention effect. The reasoning, not just the table, is the result.

---

## 1. Voxel sweep / Table 3 — PREDICTION

**Table 3 residuals reproduce to within the seed spread (0.80 % mean / 1.98 %
worst). Any tissue moving further is the finding.**

Basis — not assumption, measurement. Sim-1's three convention-fixed tubes are
the direct analogue: the same perturbation, the same estimator family, measured
on three independent materials.

| tube | Δ yield | Δ tau |
|---|---|---|
| water_shibuya2020_n2 | +32.7 % | **−0.05 %** |
| hepatic_avachat2024 | +32.2 % | **−0.15 %** |
| muscle_avachat2024 | +31.8 % | **+0.26 %** |

The estimator separated tau from I3 under exactly this perturbation, three
times.

**What a failure would mean.** If the voxel-level truncation inversion assumes a
fixed I3 anywhere, this is where it surfaces — as tissue residuals that track
the +31 % yield shift instead of staying flat. That would matter more than the
number it changes, because it would mean the estimator is not I3-invariant and
every derived tau in the paper inherits an assumption nobody declared.

## 2. Migration budget — PREDICTION

**Unaffected.** The 22.6 % within-1 mm figure and the per-material net fluxes
are truth-level geometry: where positrons stop and annihilate, which is set by
transport and the material map, not by the probability that the stopped positron
forms positronium. I3 enters only *after* the annihilation site is fixed.

**To be confirmed by measurement, not by this argument** — one number checked
against the May value.

## 3. Blind chain and unmixing — PREDICTION, lower confidence

The blind products are **stale until re-derived**, and every number in that
subsection is provisional in the same sense as Table 3.

| product | prediction | confidence |
|---|---|---|
| activity reconstruction | unchanged (activity is I3-independent) | high |
| per-voxel tau maps | unchanged, per §1 | high |
| blind region lifetimes 1.80-2.45 ns | unchanged | medium |
| **K = 2 selection** | **may change** | low |
| **E0 = 1.65, E1 = 2.9 ns** | **may change** | low |
| AUCs 0.997 / 0.944 / 0.956 | survive | medium-high |

The unmixing operates on **delay spectra**, and a +31 % change in the o-Ps
weight changes the amplitude of the o-Ps component relative to the prompt and
para terms. That can move component amplitudes, the model-order selection, and
the fitted component lifetimes — even where per-voxel tau is invariant. K
selection is a split-half reproducibility threshold, so an amplitude change can
tip it either way.

The AUCs measure **spatial contrast** between abundance maps rather than
absolute amplitude, so they should survive. That is a prediction, not an
assumption, and will be reported as measured.
