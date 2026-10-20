# P5 — Quadra energy-overflow prompt selection. Registered BEFORE the run.

Written 2026-08-24, after the referee report, before any overflow-window
extraction exists.

## The fidelity gap being tested

The Biograph Vision Quadra digitises photon energy only up to 726 keV. Every
photon above that lands in one overflow bin, so the instrument CANNOT
distinguish a 1157 keV Sc-44 prompt from anything else above 726 keV. The
published experiment (Mercolli et al., arXiv:2506.13460) therefore selects
prompts in a 720-735 keV bin that is inclusive of everything above 726 keV.

Our simulation instead applies [1100, 1200] keV to the digitised (10.1 % FWHM
blurred) energy. That is an energy discrimination the modelled scanner does
not possess.

## What is already measured (singles level, C1 NEMA at 388 ps, 178,490,874 singles)

    [1100,1200] keV prompt candidates    34,298,435   (19.216 % of singles)
    >= 726 keV   prompt candidates       38,672,795   (21.667 % of singles)
    ratio                                1.128x

    newly admitted 726-1100 keV           1,694,962   (0.950 %)
    newly admitted 1200-1300 keV          2,679,456   (1.501 %)

Note the 1200-1300 keV band: a 1157 keV photon cannot deposit more than
1157 keV, so those events are the blur tail. At 10.1 % FWHM referenced to
511 keV the resolution at 1157 keV is about 6.7 % FWHM (sigma ~ 33 keV), so
[1100,1200] is only about +/-1.5 sigma about the peak and is ALREADY cutting
genuine prompts.

## Predictions

P5-a  The genuine-triple fraction (all three singles sharing an eventID)
      FALLS under the overflow window. Predicted fall: 5 to 15 percentage
      points. Basis: the added candidates are 1.13x in number and are
      dominated by Compton continuum and blur tail, which are uncorrelated
      with the annihilation pair more often than a full-energy prompt is.

P5-b  The recovered truncated-mean tau FALLS, because a higher random
      fraction phi pulls tau_obs toward the pedestal by the registered
      pedestal law tau_obs = [m1 - phi*m_ped]/(1-phi) with m_ped = 2.0 ns.
      For tau < m_ped the pull is downward. C1 closest is 1.1202 ns.
      Predicted new value: 1.03 to 1.12 ns.

P5-c  sigma(tau) does NOT improve despite ~13 % more candidates, because the
      added events are disproportionately background. Predicted: sigma within
      a factor 0.9 to 1.3 of the C1 value (0.0007 ns paired).

## The falsifier that matters most

If the genuine fraction barely moves and tau barely moves, then our digitizer
CANNOT reproduce the background the real instrument reports (pSBR 12.6,
>20 % background inside spheres). That would mean the window is NOT the
dominant fidelity gap, and the missing physics is PILE-UP and DEAD TIME,
which the production digitizer does not model at all (no pileup, no deadtime
module in any production macro).

That outcome is the more likely one and it is the more important result:
re-windowing alone would then be insufficient, and the manuscript could not
claim a Quadra-faithful prompt selection merely by changing the window.

## Method

The frozen T6 tool (nema_extract.py) is NOT modified. Instead a derived
singles file is written with in_prompt recomputed from the stored
energy_keV; every other column is copied byte-for-byte. The unmodified
extractor is then run on both files. Control: re-running the derived-file
path with the ORIGINAL [1100,1200] window must reproduce C1 closest =
1.1201616693843575 exactly, or the derived file is wrong and nothing
downstream is trusted.
