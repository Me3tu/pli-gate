# V2 — registered prediction, written BEFORE the code change and the run

## What changes

`G4PsAtRestProcess::SamplePsLifetime` and `G4PositroniumAnnihilation::SamplePsLifetime`
both currently do:

    if (G4UniformRand() > prob) { channel = kDirect; return 0.0; }

i.e. direct annihilation is emitted at **exactly zero delay**. The new arm draws
it from an exponential with tau_direct = 388 ps (mercolli2025iodine's fit
constant; their Table 3 fixes tau_2 = 0.388 ns for direct annihilation).
p-Ps stays at 125 ps. Switch: `PLI_DIRECT_TAU_PS`, default **0 = current
behaviour, bit-exact** (the extra random draw is short-circuited, so the
baseline stream is unchanged).

## Why this could matter

The matching artefact is decided before the 1 ns censor, so a change to the
early-time model is not obviously harmless. Concretely, direct annihilation is
the majority channel, and at 388 ps a non-trivial tail leaks past the censor:

    P(delay > 1 ns | tau = 388 ps) = exp(-1/0.388) = 7.6 %

against zero leakage today. Those are **genuine** (same-eventID) triples, so
they enter `n_true_sig`, the denominator of phi.

## Predictions

1. **Absolute bias worsens under both rules.** The leaked direct events are
   short-lived relative to o-Ps, so they steepen the fitted tail and pull
   tau-hat further below the configured value. Applies to closest AND
   all-prompts.

2. **The artefact (the difference between rules) SHRINKS, and does not vanish.**
   phi = over/n_oPs; the numerator is set by the matching rule and should be
   unchanged, while the denominator gains the leaked direct events. By the
   pedestal law Delta|bias| = 1.15 * Delta_phi, so the artefact should fall by
   roughly the fractional growth of the fit-window genuine count.

3. **Quantitatively**: I expect the artefact to retain **more than half** its
   current magnitude. If it drops below half, or grows, that is a finding about
   the early-time model and outranks the rest of the block.

## The statistic

A 2x2 of extracted tau: {closest, all-prompts} x {direct at 0, direct at 388 ps},
one randoms-bearing configuration (water cylinder, 7 MBq), production
statistics. The quantity of record is the **difference of differences**

    DiD = [tau_closest - tau_all]_(388ps)  -  [tau_closest - tau_all]_(0ps)

with its uncertainty. DiD = 0 means the artefact is insensitive to the
early-time model; the sign of a non-zero DiD says which way.

## Failure modes, all reportable

* artefact unchanged within statistics -> the simplification is harmless and the
  paper says so with a number
* artefact shrinks as predicted -> quantify, and the pedestal law's denominator
  interpretation gains support
* artefact grows, or changes sign -> prediction wrong, and it outranks V3-V10

**No tuning.** tau_direct is fixed at the cited 388 ps and is not adjusted.
