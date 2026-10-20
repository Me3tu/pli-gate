# Registered predictions — N2 and N3

Written **before** `measure_oversubtraction.py` was run on any phantom.
Timestamped by the git-less convention of this project: this file precedes the
first entry in `run.log` for stage `OVERSUB`.

## N2 — the over-subtraction measurement

Definitions, fixed here and not revisited:

* fit window: `TAU_MIN_NS = 1.0 < tau < W_PROMPT_NS = 5.0`, the exact window
  `subtract_and_mle` fits
* a triple is **genuine** iff its two 511 singles and its prompt single all share
  an `eventID` within the shard. eventID labels only; it never selects.
* `n_sub` = delayed-window count in the fit window (what the estimator removes)
* `n_true_acc` = accepted signal-arm triples in the fit window that are not genuine
* `n_true_sig` = accepted signal-arm triples in the fit window that are genuine
* over-subtraction = `n_sub − n_true_acc`, error `sqrt(n_sub + n_true_acc)`

**Prediction (the authors', restated):**

1. Under **all-prompts**, `n_sub` and `n_true_acc` agree within counting error.
2. Under **closest-prompt**, `n_sub` **exceeds** `n_true_acc`.

If (2) fails — if the delayed estimate matches or undershoots the truth under
closest-prompt — the over-subtraction hypothesis is wrong and N3/N4 do not run.

## N3 — does phi explain the bias?

    phi = (n_sub − n_true_acc) / n_true_sig

**Prediction:** measured |bias| is a single monotonic function of phi, with all
four phantoms falling on **one curve**.

* one curve → mechanism established quantitatively; only the functional form
  remains open (N4)
* one curve per phantom → a second variable is in play and phi alone is
  insufficient; report which

I3 and phi are not independent — I3 sets `n_true_sig`, the denominator. **Both
correlations will be reported** so it is visible which is doing the work. A
partial correlation of |bias| with phi at fixed I3, and with I3 at fixed phi,
will be given where the sample allows.

## What would make this a failure worth reporting

Any of: n_sub < n_true_acc under closest-prompt; phi uncorrelated with bias;
phantoms on separate curves; or phi correlating only because it inherits I3.
All four are reportable outcomes and none will be repaired by adjustment.

**No tuning.** The fit window, the truth definition, the estimator and the
denominator are all fixed above, before any number was produced.

---

## N3b — the difference test (registered before the Steinberger all-prompts arms landed)

N3 found that phi collapses the closest-prompt data onto one curve, but that
the same curve **under-predicts the all-prompts bias by a median factor 2.5**.
A bias floor survives phi -> 0. Fitting a power law to the closest arm alone
absorbs that floor, which is why its exponent came out 0.656 rather than 1.

The floor is phi-independent by construction, so it cancels in a difference
taken between the two matching rules on the SAME data:

    delta_bias = |bias_closest| - |bias_all|          (both measured)
    delta_phi  = phi_closest   - phi_all              (both measured)

Every quantity is measured. There is **no fitted constant**, which is the point:
the previous test needed one and I will not introduce it.

**Prediction:** delta_bias is a single monotonic function of delta_phi, common
to all four phantoms, passing through the origin. If it is also linear, the
slope is the conversion from over-subtracted fraction to lifetime bias and is
the quantity N4 must reproduce analytically.

**Failure modes, all reportable:**
* phantoms on separate lines -> a second variable modulates the over-subtraction
  term itself, not just the floor
* non-zero intercept -> something other than over-subtraction differs between
  the rules
* no correlation -> the over-subtraction is not the causal agent, and the N3
  curve was inherited from the accidental-to-signal ratio rather than from
  over-subtraction specifically (which N3 already could not exclude:
  over-subtraction alone had Spearman -0.303, p = 0.31)

Pairing rule: both arms of every difference come from the same singles file and
the same extractor invocation, differing ONLY in --matching. Sim-1's pair is
`sim1_{closest,all}_tau` at sigma_T = 90.9 ps; the Steinberger pairs are
`{disks,water}_{closest,all}_o4`; NEMA's is `out_{closest,all}` from T6.
