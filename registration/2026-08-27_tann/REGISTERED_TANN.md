# Registered before the t_ann comparison ran — 2026-08-27

## What is compared

Two reconstruction chains on IDENTICAL selected triples:

    explicit  tau = (t1 - d1/c) - (tp - dp/c)          per_event_tau, PRODUCTION
    simplified tau = 0.5*(t1 + t2) - tp                the t_ann = (t1+t2)/2 chain

Only `per_event_tau` is patched, and it is called AFTER selection, so both runs
select exactly the same triples with the same seed. Nothing else changes.

## The algebra, written down before the measurement

For a back-to-back pair emitted at p_ann, t1 = t_ann + d1/c and t2 = t_ann + d2/c, so

    tau_simplified - tau_explicit = (d1 + d2)/(2c) - dp/c  =  -Delta_geo

which is exactly the path-length term whose spread the acceptance kernel already
carries (sigma_geo 305-353 ps depending on geometry; see sigma_geo_sim1.json).
So the difference is NOT a bias in the lifetime; it is a per-event offset with
zero delay dependence.

## Registered expectation (the authors', restated)

Consistent within statistics, because the path spread is delay-independent and
enters only the window-edge acceptance that the estimator models.

    P-T1  |delta tau_hat| < 2 x its paired sigma on every VOI.
    P-T2  The per-triple difference has a spread of order sigma_geo (300-360 ps)
          and a mean of order +10 to +40 ps, matching Delta_geo's own statistics.
    P-T3  Any residual delta tau_hat does NOT correlate with tau_cfg, since the
          mechanism is delay-independent.

If P-T1 fails the acceptance treatment needs revisiting, and that is reported as
a finding rather than absorbed.
