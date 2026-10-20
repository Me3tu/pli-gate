# G6 — channel composition of the fit window. Registered before any tagging.

## Why this is now the only live route

The estimator route to the 9-16 ps chain offset is closed: pairing the two
estimators shrinks sigma only to 0.72-0.81 because their difference lives in
the tail above the 5 ns cut. Additive delays are excluded analytically. What
remains is a change in the POPULATION the estimator reads.

## The test

Every annihilation is emitted through exactly one channel: direct (0 ps by
default), p-Ps (125 ps), o-Ps 2-gamma or o-Ps 3-gamma (both at tau_eff). The
physics writes the sampled lifetime per annihilation, so the fit-window
population can be decomposed by channel without any modelling.

Measure, per region and separately on axis and off axis, the fraction of the
population above the 1 ns censor contributed by each channel.

## The arithmetic that makes this a real test

Above the censor an exponential is memoryless, so the estimator returns the
count-weighted mean of the component lifetimes. A contaminant of fraction f
with lifetime tau_c shifts the answer by

    delta = f * (tau_c - tau_oPs)

For a SHORT contaminant (p-Ps at 0.125 ns, tau_oPs 1.84 ns) the shift is
NEGATIVE: f = 1 % gives -17 ps. The observed offset is POSITIVE, +9 to +16 ps.
A short-lived contaminant therefore CANNOT produce it. To get +12 ps needs
either a LONG component, tau_c > tau_oPs, at about 1 % if tau_c were ~3 ns, or
a deficit of short-lived events.

**Registered prediction.** I expect the fit-window population to be
overwhelmingly o-Ps with p-Ps at the 1e-4 level, because p-Ps leakage past a
1 ns censor is exp(-8) = 3.4e-4 of its population. If so, no channel
contaminant of the required size exists and the channel hypothesis is
FALSIFIED, leaving the cause outside the emitted population entirely. If
instead an unexpected long-lived component appears at the percent level, that
is the cause and it is directly visible.

**Secondary prediction.** Composition will be the same on axis and off axis to
within counting error. If it is not, the on/off-axis difference in the offset
(+9.25 vs +15.62 ps) has a composition explanation.

**No tuning.** Channel labels come from the physics, not from a fit.
