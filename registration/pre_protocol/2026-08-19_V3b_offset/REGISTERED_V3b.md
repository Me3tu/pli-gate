# V3b — registered before the run: where does the +15.6 ps offset enter?

## The observation this tests

V3 measured a common offset of **+15.62 +- 1.66 ps** between the recovered
oracle lifetime and the configured one, flat in pO2 across 0-3000 mmHg, and not
radial. V2's arm A shows **+8.54 +- 2.8 ps** in its best-measured VOI, on axis.
It is 0.85 % at tau = 1.84 ns and larger than the whole physiological oxygen
signal, so it matters more than the thing V3 set out to measure.

## The discriminator

The chain has two halves and the offset is in one of them:

    sampler:        SamplePsLifetime draws tau ~ Exp(tau_eff)
    reconstruction: detector hit times -> prompt match -> delay -> estimator

`PLI_LIFETIME_LOG=1` writes one record per annihilation,
`PLI_LIFETIME:<tau_sampled>:<material>:<x>:<y>:<z>:<global_time>`, where field 1
is the **lifetime the sampler emitted** — before positron transport, before
photon time of flight, before the 90.9 ps timing smear, before prompt matching.

Applying the production estimator's own arithmetic, `mean(tau | tau > 1 ns) -
1 ns`, to field 1 and to the reconstructed delays from the SAME run separates
the two halves cleanly. Nothing else differs between the two numbers.

## Predictions

**I expect the sampler to be exact and the offset to be added downstream.**
`SamplePsLifetime` draws `-tau * log(U)`, which is exact by construction, so
there is no room in it for a 0.85 % bias. Specifically:

1. **Sampled**: `mean(tau_sampled | > 1 ns) - 1 ns` = 1.840 ns, within its
   statistical error. If it is not, the sampler is emitting the wrong
   distribution and that outranks everything else in this block.
2. **Reconstructed**: the oracle on the same run returns ~1.856 ns, reproducing
   the offset.
3. Therefore the offset is **added by the reconstruction**, and the candidates
   are positron thermalisation and transport time before annihilation, photon
   time of flight, and the prompt reference.

## What each outcome would mean

* sampled exact, reconstructed offset -> downstream. Next step is to split
  transport from TOF, which the record's global-time field supports.
* both offset -> the sampler. A 0.85 % bias in an exponential draw would be a
  serious bug and would invalidate configured-vs-recovered comparisons
  throughout the project.
* neither offset -> the offset is specific to the V3/V2 configurations rather
  than general, and the comparison has to be repeated on those exact runs.

**No tuning.** Nothing is adjusted to make either number land on 1.84.
