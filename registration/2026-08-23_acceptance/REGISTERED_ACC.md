# Acceptance term + D2 unification — registered 2026-08-23 16:50, before any fit

Both tests below were specified by the authors in the 2026-08-23 decision. The bars are
transcribed here verbatim so they cannot drift, and the derivations that follow
were done before any acceptance-corrected fit was run.

## Part 1 — the acceptance function A(tau), zero free parameters

### Route 1, analytic

Derivation, done here so the constants are visibly not fitted. Write the box the
extractor applies, |t_prompt - t_pair_mean| <= W with W = 5 ns, in terms of the
quantity that is actually histogrammed.

With true emission times and per-single Gaussian jitter eps_i:

    t_i = t_ann + d_i/c + eps_i        (i = 1, 2, the two 511s)
    t_p = t_pe  + d_p/c + eps_p        (the prompt)

The extractor forms tau = (t_1 - d_1/c) - (t_p - d_p/c) and the box variable
u = t_p - (t_1 + t_2)/2. Substituting and cancelling:

    tau_meas = (t_ann - t_pe) + eps_1 - eps_p
    u        = -tau_meas + Delta_geo + (eps_1 - eps_2)/2

where Delta_geo = [d_p - (d_1 + d_2)/2]/c is purely geometric. So the box
|u| <= W becomes, in the histogrammed variable,

    |tau - S| <= W,     S = Delta_geo + (eps_1 - eps_2)/2

and therefore

    A(tau) = Phi((W - tau)/sigma_e) * Phi((tau + W)/sigma_e)

with the window edges at exactly +/-5 ns BY CONSTRUCTION, and

    sigma_e = sd[(eps_1 - eps_2)/2] = sigma_t / sqrt(2) = 90.9 / 1.4142
            = 64.3 ps                                    (+ a small Delta_geo term)

sigma_t = 90.9 ps is the production per-single timing sigma, fixed long before
this test. **No constant in A(tau) is fitted.**

### Route 2, empirical

The delayed/randoms arm is flat by truth, so its measured shape IS the
acceptance. Measured in A1 on 697,569 pooled counts: flat to +/-1.2 % across
[-4.5, +4.5] ns, then 0.712 at the low edge bin and 0.766 at the high edge bin.
Normalise by the flat interior and use directly.

### Registered agreement test

R1  Routes 1 and 2 must agree. Concretely: fitting sigma_e to route 2's measured
    shape must return a value consistent with the analytic 64.3 ps. That fitted
    sigma_e is used ONLY as the agreement statistic; the re-fit uses the analytic
    value, so the zero-free-parameter claim holds.

## Part 2 — the gate, re-run THROUGH the acceptance

A gate that cannot represent the failure mode it guards against is not a gate.
Synthetic spectra are regenerated WITH A(tau) applied, then fitted by the
acceptance-corrected model. Same four conditions as B2:
1 % of truth and 3 sigma at 1e6, o-Ps share within 0.02, 3 sigma at 1e5.

## Part 3 — the re-fit, bars unchanged

R2  tau_oPs on the 388 ps water arm lands in **[1.75, 1.93] ns**, both rules.
R3  The endpoint sweep that swung **407-625 ps** collapses to **under 3x the
    statistical error** (statistical error is 9-12 ps, so under ~30-36 ps).
R4  NOTHING is retuned per sample. If it lands only with per-sample adjustment,
    STOP and book the failure.

## Part 4 — the D2 unification test

D2 (V3b, `findings_V3b.md`): the tube residual falls with configured lifetime,

    slope        = -19.11 +/- 6.27 ps/ns   (3 sigma from flat)
    zero crossing near tau = 2.01 ns
    +8.1 ps at 1.77 ns  ->  -14.6 ps at 2.72 ns

The proposed mechanism: an acceptance deficit at the 5 ns edge removes proportionally
MORE tail from LONGER lifetimes, so tau_hat is pulled down more at large tau,
giving a residual that FALLS with tau.

Propagate the measured A(tau) through the production truncated-mean estimator as
a function of tau, with no new constants, and compare.

R5  **The SLOPE is the zero-parameter prediction.** Predicted d(resid)/d(tau)
    must match -19.11 +/- 6.27 ps/ns within the error.
R6  **The zero crossing is NOT zero-parameter and is recorded as such.** The
    acceptance mechanism alone loses tail at every tau, so it predicts a residual
    that is negative everywhere. A crossing at 2.01 ns requires a positive
    constant added on top -- and one is already measured independently: the
    +9 to +16 ps reconstruction offset. The crossing is therefore a CONSISTENCY
    check on that known constant, not a free fit, and it will be reported that
    way.
R7  Falsifier: if the predicted slope is wrong by more than the 6.27 ps/ns error,
    the edge deficit is a mechanism but not THIS one. Both facts go on the page.

R8  Unchanged regardless of outcome: the deficit does NOT explain the pull
    outliers (T8 myxoma has the larger deficit and one ninth the residual), so
    the pulls page stands as written.
