# B2 registered before the bootstrap runs — 2026-08-26

## What B2 measures

The paired uncertainty on the rule-difference DiD = |bias_closest| - |bias_all|
for every region of the C1 (NEMA IQ) and C2 (Sim-1) 388 ps arms, and for the
corresponding frozen 0 ps arms, by resampling SHARDS coherently: one draw picks
a shard set and BOTH rules are recomputed from exactly those shards.

Definitions are taken verbatim from the frozen tools, not restated:
  bias    = 100*(tau_trunc - tau_cfg)/tau_cfg          test_delta.py:56-57
  DiD     = |bias_closest| - |bias_all|                test_delta.py:62
  tau     = truncation_correct(mean(sig-rnd) - TAU_MIN) v4d_paired.tau_from_hists

The production number this replaces is did_sigma = hypot(sigma_c, sigma_a),
c2_result.py:67 — the INDEPENDENT combination, which assumes rho = 0.

## Prior

V4(d) measured, on the five Steinberger regions:
    rho     +0.967 to +0.981   (median +0.972)
    ratio    0.141 to  0.187   (median  0.176)  paired / independent

## Registered predictions

P-B2-1  rho > 0 in every region of both phantoms and both arms. A negative
        correlation anywhere is a STOP: it would mean the two rules anti-track
        across shards, which no mechanism in the pipeline explains.

P-B2-2  The paired/independent ratio is BELOW 1 everywhere, so every published
        error bar on a rule-difference shrinks. Central expectation is the V4(d)
        value 0.176; the registered band is 0.10 to 0.45. Outside that band is
        reported as measured, NOT rescaled onto it.

P-B2-3  The 388 ps arms have LOWER rho than the 0 ps arms. Reason: at 388 ps the
        direct component is displaced into the fitted window, so the two rules
        see partially different event populations rather than the same one
        re-matched. If rho instead rises, the displacement is not doing what the
        early-time model says it does, and that is a finding.

P-B2-4  Applying the paired sigma to B3 raises chi2/dof of the pedestal law
        ABOVE 1 and so yields a NON-ZERO sigma_sys. B3 with independent errors
        gave chi2/dof = 0.221; shrinking sigma by the V4(d) factor 0.176 would
        multiply chi2 by ~32, which massively overshoots 1. So the honest
        prediction is: a real systematic floor appears, and it is the DOMINANT
        term in the residual. B3 must then be re-run and re-reported.

## What is NOT allowed

No rescaling of any sigma onto the registered band. No dropping of a region
because its rho is inconvenient. If P-B2-1 fails anywhere, surface and STOP.
