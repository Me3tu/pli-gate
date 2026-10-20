# Observation registered, NOT chased: the Mercolli oracle sits ~300 ps high

Found 2026-08-23 while preparing B3. Recorded now, before any test, so that a
later explanation is checked against a prediction rather than fitted to the
number. Decision 2 of the 2026-08-23 evening block says second-order
mis-specifications are registered and not chased tonight; this is registered
under that rule.

## The observation

The oracle estimator is `mean(tau | tau > 1 ns) - 1 ns`, no upper cut, on
truth-tagged eventID-grouped triples (verified at
`steinberger_extract_ap.py:258`). Above a censor an exponential is memoryless,
so the truncated mean of a mixture is the count-weighted mean of the component
lifetimes, weighted by what each component puts ABOVE the censor. That is the
same algebra V2 used, and nothing in it is fitted: I3 and tau come from the
material database, tau_direct = 388 ps from the run script, tau_pPs = 125 ps.

| tube | I3 | f_direct above censor | predicted | measured | measured - predicted |
|---|---|---|---|---|---|
| T1_dry | 0.211 | 0.278 | 1.9273 | 2.2076 | +280 ps |
| T2_wet | 0.217 | 0.275 | 1.8246 | 2.1332 | +309 ps |
| T3_gelatin | 0.222 | 0.272 | 1.7572 | 2.0563 | +299 ps |
| T4_water | 0.268 | 0.240 | 1.4763 | 1.8168 | +341 ps |

**The measured oracle is 280-341 ps ABOVE the prediction on all four tubes**,
by a roughly common amount. Less of the 388 ps direct component survives the
1 ns censor than the configured intensities imply.

## Already excluded

- **Not a wrong direct-tau setting.** The 388 ps value did take effect: this
  run's clinical tau_trunc is 1.13-1.19 ns, matching B5's 388 ps water arm
  (1.146) and nothing like its 0 ps arm (1.672).
- **Not my oversized VOI admitting air.** The oracle assigns by TRUE SOURCE
  position (`voi_masks(src, vois)` at line 214), and every source is inside a
  tube by construction, so the 15 mm VOI cannot pull in air-borne decays.
- **Not the oracle having an upper cut** that my formula omitted; it has none.

## Candidates, in the order I would test them

C1  The I3 convention actually applied by the loader. If the shipped I3 is
    consumed as total formation rather than as P(o-Ps), the direct share
    changes. NOTE this moves the prediction the WRONG way (more direct, lower
    prediction), so on its own it makes the gap worse, not better.
C2  A selection effect against direct annihilation in the genuine-triple
    requirement that I have not identified. Direct and pick-off both give a
    back-to-back 511 pair, so there is no obvious one -- which is why this needs
    a measurement rather than an argument.
C3  Positron escape from the 6 mm tubes into air, where I3 = 0 so everything is
    direct. This ADDS direct-like counts and also moves the gap the wrong way.

C1 and C3 both push the wrong direction, which is itself informative: whatever
is happening REMOVES short-lived counts from the censored sample rather than
adding long-lived ones.

## Falsifier for whatever is proposed

The gap is nearly constant in ps across four tubes spanning tau 1.82 to 2.52 ns
and I3 0.211 to 0.268. Any mechanism offered must reproduce that near-constancy,
not just the sign.

**Second falsifier, added 2026-08-24 at the authors' request — this is the harder
one.** The SAME parameter-free mixture algebra CLOSED on the V2 water dataset:
`findings_V2.md` reports predicted against measured oracle at +0.23 %, +0.07 %
and -0.34 % across three VOIs. So the algebra is not simply wrong.

Any proposed mechanism must therefore explain BOTH:
  (a) a ~300 ps uniform excess on the Mercolli tubes, and
  (b) why the identical algebra closed to better than 0.4 % on V2 water.

Candidate discriminants, none tested: a different oracle definition between the
two extractions, different censor handling, different geometry (6 mm tubes in
air against a 20 mm water cylinder), or a real difference between the runs.
**A mechanism that explains one and not the other is wrong.**

## What this does NOT affect

Nothing on the B3 or B4 pages depends on the oracle being correct: those pages
compare our CLINICAL estimator against Mercolli's published values under their
model. The oracle is reported there as a diagnostic with this gap stated.
