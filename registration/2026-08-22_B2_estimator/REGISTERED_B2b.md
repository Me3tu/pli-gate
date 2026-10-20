# B2b — the gated estimator on REAL data. Prediction registered before it runs.

Written 2026-08-22 ~19:20, after B2's synthetic gate passed and BEFORE the
estimator has seen any simulated spectrum.

## What is being tested

B5 showed that with a realistic 388 ps direct component the truncated-mean
estimator returns 1.146 ns (closest-prompt) and 1.196 ns (all-prompts) on the
Steinberger water cylinder, against a configured 1.84 ns and a published
1.815 +/- 0.013 ns. The claim on that page is that this is a property of the
SIMPLE ESTIMATOR and not of the simulation.

That claim is testable right now, on the same events, with no new simulation:
fit the same randoms-subtracted spectrum with the three-component model that
passed its gate today.

## Predictions

P1  **If the claim is right**, the three-component fit recovers tau_oPs within
    5 % of the configured 1.84 ns, i.e. inside [1.75, 1.93] ns, on the 388 ps
    arm. It must do so under BOTH matching rules.

P2  **The falsifier.** If the three-component fit also returns something near
    1.15 ns, then the deficit is not the estimator and page 19's claim is wrong.
    In that case I report the claim as refuted and do not adjust the model to
    rescue it.

P3  On the 0 ps arm the three-component fit should agree with the truncated mean
    to within its own uncertainty, because at zero delay there is nothing for the
    third component to absorb. A large disagreement there would mean the fit is
    finding structure that is not present, and would undercut P1 even if P1
    passed.

## Fixed, not fitted

tau_pPs = 125 ps and tau_direct = 388 ps are held at the values the block
specifies, exactly as in the gate. On the 0 ps arm tau_direct is held at 388 ps
as well, deliberately: the point of P3 is to see whether a component that should
have no counts to claim stays harmless. Its intensity is free and may go to zero.

## What this cannot show

Both arms come from our own simulation, so this tests the estimator against our
generator, not against nature. Mercolli's own samples are what test the model.
