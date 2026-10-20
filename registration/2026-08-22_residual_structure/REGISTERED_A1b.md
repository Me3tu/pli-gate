# Edge deficit in the delayed window — hypothesis registered before the test

Written 2026-08-22 ~18:40, after measuring the shape and BEFORE looking at the
low edge of the histogram. Nothing below may be edited after the test runs.

## The measurement that prompted it

The subtracted (delayed-window) counts are flat over [1.0, 4.5] ns to
chi2/dof = 1.33 median across the 14 populated regions, and then fall by
21.9 +/- 4.7 % in the final 0.5 ns bin [4.5, 5.0]. That deficit, not a tilt,
is what produces the -62 ps first-moment shift reported on page 8: a uniform
density with a 22 % hole in its last 12.5 % has first moment 1.951 ns above the
censor against 2.000 for a true uniform, and the measured value is 1.938.

## Candidate mechanism

`clinical_per_shard` selects the delayed prompt with a hard box on RAW detector
time, `|t_prompt - (t_pair_mean + DELAY_NS)| <= W_PROMPT_NS`, but the quantity
histogrammed is `per_event_tau`, which subtracts a per-event time-of-flight
correction `d/c` from each arm. The box is therefore sharp in one variable and
blurred in the other by the spread of that correction, so acceptance falls off
near the window edge rather than cutting cleanly.

## Prediction registered NOW, before the low edge is examined

P1  If this is the mechanism, the deficit is a property of the BOX EDGE and not
    of the long-delay side, so it must appear at BOTH edges. The stored
    histogram spans [-5, 8] ns, so the low edge near tau = -5 ns is on disk and
    has never been looked at. The delayed-window density must fall away near
    -5 ns by a comparable fraction.

P2  A deficit at the high edge ONLY refutes the mechanism. In that case the
    cause is something that acts on long delays specifically, and I will say so
    and register a new candidate rather than reinterpret this one.

P3  The size of the edge roll-off is set by the spread of the time-of-flight
    term, which is a geometric property of the region, so it must NOT correlate
    with the region's residual in the identity test. If it does correlate, the
    mechanism is entangled with the outliers and the claim on page 8 that
    hypothesis (3) fails would need revisiting.

## What is NOT being done

No code is changed. The frozen numbers were produced by this selection and the
paper must describe what it did. This is a measurement of a known-shape
assumption, not a repair.
