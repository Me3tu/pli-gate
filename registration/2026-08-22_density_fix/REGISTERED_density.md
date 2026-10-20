# Density correction — registered before the confirmation run

## What changed

`XAD4_Dry` 1.02 -> 0.67 g/cm3 and `XAD4_Wet` 1.08 -> 0.67 g/cm3, to the
measured rho_tot of Table 1. Dry also gains the 0.4 ml of [124I]NaI solution
every tube contains (C 0.8103 / H 0.0924 / O 0.0973 by mass, from 3.25 g
polymer + 0.40 g solution). Wet's composition is unchanged. Gelatine and water
already matched.

## Predictions

1. **tau is untouched.** The o-Ps lifetime is read from the materials JSON and
   is not a function of density anywhere in the code. The extracted tau must
   agree with the pre-fix run to within counting error. This is verified, not
   assumed, because "by construction" claims are exactly the ones that turn out
   to have an unnoticed coupling.
2. **Scatter increases, and the singles yield falls.** Less material means
   fewer interactions in the sample: fewer annihilations occur (the source is
   confined to the tube volume, so decay count is unchanged, but attenuation
   and in-sample scatter drop). The scattered fraction reaching the detector
   from within the tube should FALL, and the unattenuated 511 keV flux should
   RISE. Direction is predicted; magnitude is not.
3. If tau moves by more than counting error, the density is coupled to the
   lifetime somewhere it should not be, and that outranks the exposure run.

**No tuning.** Densities are Table 1's measured values, not adjusted to
improve anything.
