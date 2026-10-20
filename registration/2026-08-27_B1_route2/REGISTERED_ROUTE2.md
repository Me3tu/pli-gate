# Route 2 registered before it ran — 2026-08-27

## The rules I am bound by (from DERIVATION_lever_2026-08-26.md)

> Do NOT implement the closed form. Verify the DEFINITION directly, with
> machinery of your own choosing [...] Two evaluation points minimum
> (tau_hat ~= 1.14 and ~= 1.68). Agreement bar: 1 % on L at both points.
> If your 0.1996 came from a different evaluation point or an uninverted
> mean, this procedure will show exactly where it diverged — report the
> divergence point, not just the verdict.

So: the expression (m_ped - m1(tau_hat)) / m1'(tau_hat) is NOT written anywhere
in my code. Neither is m1' in any form. What I use is:

  1. a window population S on [c, C] = [1, 5] ns, built from its density;
  2. a flat component mixed in at fraction phi, per the contamination model
     m_meas(phi) = (1 - phi) m(S) + phi m_ped;
  3. **the production `truncation_correct` imported from
     steinberger_extract_ap.py**, not a reimplementation;
  4. L = d tau_hat / d phi at phi -> 0 by central difference, with the step
     refined until the estimate stops moving.

Two independent constructions of the population, which must agree:
  * deterministic quadrature of the density (noise-free), and
  * Monte-Carlo sampling of the same density (noise, but no quadrature
    assumptions).

## Targets, from route 1 (the authors', already done)

    L(1.1366) = 1.6146      pure exponential
    L(1.6842) = 2.0038      pure exponential
    L(mixture 1.80 + 26 % x 0.388, tau_hat = 1.2374) = 1.6757
    ratio L_388 / L_0 = 0.8055

## Registered predictions

P-R2-1  Both evaluation points agree with route 1 to better than 1 %.
P-R2-2  The mixture row agrees too. This is the row that tests the identity
        m(S) = m1(tau_hat0) rather than just the algebra, so if anything
        breaks it should break here.
P-R2-3  The 0.1996 I previously quoted is recoverable as a SPECIFIC wrong
        step, and I will name which. My leading hypothesis, written down
        before testing: I evaluated the lever at a CONFIGURED lifetime
        (1.80 ns) instead of the arm's OBSERVED tau_hat, and/or used the raw
        window mean without inverting it. If neither reproduces 0.1996 I say
        so rather than inventing a story.

If route 2 disagrees with route 1 by more than 1 %, C2-b stays UNSCORED and I
report the disagreement. Only if it agrees do I go on to score C2-b.
