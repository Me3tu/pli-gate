# A6 — Moskal arXiv:2605.00021 against our oxygen response. Prediction registered first.

Written 2026-08-22 ~18:45, after reading his Table II and BEFORE running our own
tao_eldrup module at the same pO2 points.

## His numbers, transcribed verbatim from Table II

Definitions are his: Delta tau_oPs = tau(pO2 = 10 mmHg) - tau(pO2 = 50 mmHg),
i.e. hypoxia minus physoxia. Delta lambda_O2 = lambda(50) - lambda(10).

| substance | Delta lambda_O2 | Delta tau_oPs | Delta R_oPs-3g/2g | R_oPs-3g/2g at pO2=0 |
|---|---|---|---|---|
| Water | 1.4 us^-1 | 5 ps | 0.000045 | 0.0164 |
| Isopropanol | 14.0 us^-1 | 182 ps | 0.0015 | 0.0298 |
| Cyclohexane | 20.7 us^-1 | 187 ps | 0.0015 | 0.0254 |
| Isooctane | 35.6 us^-1 | 444 ps | 0.0036 | 0.0309 |
| Adipose | 7.0 us^-1 | 48 ps | 0.0004 | 0.0218 |

His oxygenation definitions: physoxia = median tissue 30-57 mmHg, mild hypoxia
< 20, moderate < 10, severe < 1 mmHg. He evaluates at 10 vs 50 mmHg.

## The structural difference, stated before any number is computed

Our model uses ONE oxygen quenching constant for every material,
kappa_O2 = 2.45e-5 /ns/mmHg, so our quenching rate increment over a given pO2
step is the SAME for all materials and our Delta tau depends on the material
only through tau^2 (since Delta tau ~= tau^2 * Delta lambda for small steps).

His Delta lambda_O2 varies from 1.4 to 35.6 us^-1 across his five substances, a
factor of 25, because he carries a material-specific oxygen solubility.

## Predictions registered NOW

P1  For WATER the two should roughly agree. Our Delta lambda over the 40 mmHg
    step is 40 * 2.45e-5 /ns = 0.98 us^-1 against his 1.4 us^-1, so our Delta tau
    should come out LOWER than his 5 ps by roughly the same ratio, order 3-4 ps.

P2  For ADIPOSE we should disagree badly and in a known direction. His
    Delta lambda_O2 is 7.0 us^-1, about 7x our material-independent 0.98 us^-1,
    so our Delta tau must be several times SMALLER than his 48 ps.

P3  The ratio (his Delta tau)/(ours) must track (his Delta lambda_O2)/(0.98 us^-1)
    to within the tau^2 factor. If it does not, the disagreement is not solubility
    and I will say so rather than assert a cause.

## What this comparison is NOT

It is not a reconciliation. Nothing in our database or code is changed to match
him. Disagreements are flagged and left standing, as specified.
