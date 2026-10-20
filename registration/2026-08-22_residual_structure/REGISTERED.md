# Residual structure of the pedestal law — hypotheses registered BEFORE any plot

Written 2026-08-22. Nothing in this file may be edited after the tests are run;
outcomes go in `findings.md`.

## The thing being explained

G3(b) replaced the independent-bootstrap error on the pedestal-law identity test
with the paired one. The point estimates did not move; the errors shrank by a
median factor of ~4.5. The consequence is that the scatter about the identity
line, RMSE = 3.74 % of tau, is no longer absorbed by the error bars: pulls now
run from -8.19 to +4.68 and 3 of 20 regions sit beyond 3 sigma.

That scatter is REAL STRUCTURE, not noise. This page asks what it is a function
of. It is not an attempt to make the law fit better: the law has no free
parameter and none will be added.

Residual is defined throughout as

    residual_i = measured_i - predicted_i        [ % of tau, DiD units ]

exactly the numerator of the G3(b) pull. Predicted comes from
`derive_pedestal_bias.predicted_bias(tau_i, d_phi_i)`, unchanged.

## H1 — the residual is a function of the configured lifetime

The lever `predicted_bias` is a nonlinear function of tau. If its tau-dependence
is imperfectly captured, the residual inherits a tau-dependence.

**Prediction registered by the authors:** "pulls sort by lifetime, with hepatic and
polycarbonate at opposite extremes of tau".

**Checked before running, as specified:**

| | region | tau_cfg | pull |
|---|---|---|---|
| most negative | T4_hepatic_avachat2024 | 2.040 ns | -8.19 |
| most positive | polycarbonate | 2.070 ns | +4.68 |

They are 0.030 ns apart — 2.6 % of the 1.56–2.72 ns range spanned by the 20
regions. They are not at opposite extremes of tau; they are on top of each
other. Pearson r(tau, pull) over all 20 regions = -0.060.

**H1 is falsified before the test is run.** The panel is still drawn, because
the reader is owed the evidence for the falsification and not only the claim.

## H2 — the residual tracks the randoms fraction

The law is exact in phi to all orders as written, but phi itself is measured per
ROI. If phi is mis-measured, or if the removed content is not exactly the
phi*N assumed, the residual grows with phi.

**Prediction registered now, before the panel exists:** if H2 explains the
scatter, |residual| increases monotonically with d_phi, and the two 3-sigma
outliers are among the largest-d_phi regions. A residual that is flat in d_phi,
or whose outliers sit at small d_phi, refutes it. No constant is fitted; the
ordering is the whole test.

## H3 — the pedestal is not flat

This is the law's ONE structural assumption, stated in the derivation docstring
as "the excess is flat in delay — accidentals are uniform in the window by
construction" and entering the arithmetic as the single line

    m_ped = WIDTH / 2.0        # = 2.0 ns

If a region's delayed-window content is not uniform over [1.0, 5.0] ns, its true
first moment m_ped* differs from 2.0 and the law is being evaluated at the wrong
pedestal.

**This is measurable with no new constant.** m_ped* is the first moment of the
delayed-window randoms histogram already stored per region in the .npz files.
Substituting the measured m_ped* for 2.0 in `predicted_bias` gives a corrected
prediction that is still parameter-free — every input is measured, nothing is
fitted.

**Predictions registered now, before the panel exists:**

- P3a  If H3 explains the scatter, the residual correlates with (m_ped* - 2.0)
       with the sign the law dictates: a pedestal whose mass sits at LONGER
       delay than uniform over-subtracts from the long side, so the recovered
       tau falls further and the measured DiD bias exceeds the prediction —
       residual POSITIVE for m_ped* > 2.0.
- P3b  Substituting m_ped* reduces the RMSE about the identity line below
       3.74 %. If the RMSE does not fall, H3 does not explain the scatter,
       and I will say so.
- P3c  The falsifier: the two 3-sigma outliers must be the two regions with the
       largest |m_ped* - 2.0|. If a region sits at 3 sigma with a pedestal flat
       to within the measurement error, H3 cannot be the explanation for it.

## Outcome rules

- All three panels are drawn regardless of outcome, including H1's.
- Nothing here is refitted. The corrected prediction of H3 uses only measured
  quantities.
- If all three fail, the page reports three failures and the residual stays
  unexplained. That is a publishable state; a fitted fudge factor is not.
