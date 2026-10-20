#!/usr/bin/env python3
"""sigma(Delta_geo) for Sim-1's own source geometry -- CONDITION 1 of the
three-component work.

Same method and the same imported detector model as
2026-08-23_acceptance/sigma_geo.py and 2026-08-24_acceptance_nema. Only the
source volume changes: Sim-1's activity is confined to eight tubes,
r = 15 mm, halfz = 25 mm, centres on a ring of radius 80 mm at 45 degree
intervals, all at z = 0 (sim1_sc44_388_timedriven.mac:156-170 and siblings).

The eight tubes are geometrically equivalent by symmetry, so one tube fixes
sigma_geo for all. That symmetry is VERIFIED here, not assumed.

REGISTERED before running (2026-08-26): Sim-1's tubes sit at radius 80 mm and
are axially short (halfz 25 mm) compared with NEMA's body (halfz 90 mm). The
axial term dominates the spread, so sigma_geo(Sim-1) should be SMALLER than
NEMA's 351.4 ps, and closer to the Mercolli value. If it comes out larger,
the geometry argument is wrong and that is a finding.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.environ.get("PLI_RES") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                 "2026-08-23_acceptance"))
from sigma_geo import hit_distance, C  # noqa -- the frozen detector model

SIG_T = 0.0909                       # ns, production single-photon timing sigma
SIGMA_T_PAIR = SIG_T / np.sqrt(2.0)
RING_R, TUBE_R, TUBE_HZ = 80.0, 15.0, 25.0
ANGLES = np.arange(8) * (np.pi / 4.0)


def sample_tube(k, n=2_000_000, seed=20260826):
    """Uniform inside tube k of the Sim-1 ring."""
    rng = np.random.default_rng(seed + k)
    cx, cy = RING_R * np.cos(ANGLES[k]), RING_R * np.sin(ANGLES[k])
    r = TUBE_R * np.sqrt(rng.random(n)); th = rng.uniform(0, 2 * np.pi, n)
    p = np.c_[cx + r * np.cos(th), cy + r * np.sin(th),
              rng.uniform(-TUBE_HZ, TUBE_HZ, n)]

    def iso(m):
        c = rng.uniform(-1, 1, m); s = np.sqrt(1 - c * c)
        ph = rng.uniform(0, 2 * np.pi, m)
        return np.c_[s * np.cos(ph), s * np.sin(ph), c]

    u = iso(n)
    d1, ok1 = hit_distance(p, u)
    d2, ok2 = hit_distance(p, -u)
    v = iso(n)
    dp, okp = hit_distance(p, v)
    ok = ok1 & ok2 & okp
    return (dp[ok] - 0.5 * (d1[ok] + d2[ok])) / C


out = {"tubes": {}}
sg = []
for k in range(8):
    d = sample_tube(k)
    s = float(d.std(ddof=1))
    sg.append(s)
    out["tubes"][f"tube_{k}"] = dict(n=int(len(d)), sigma_geo_ps=1000 * s,
                                     mean_delta_ps=1000 * float(d.mean()))
    print(f"  tube {k} (theta = {np.degrees(ANGLES[k]):5.1f} deg)  "
          f"sigma_geo = {1000*s:6.1f} ps   n = {len(d):,}")
sg = np.array(sg)
spread = 1000 * (sg.max() - sg.min())
print(f"\n  SYMMETRY CHECK: spread across the eight tubes = {spread:.2f} ps "
      f"({100*(sg.max()-sg.min())/sg.mean():.3f} % of the mean)")
if spread > 2.0:
    raise SystemExit(f"ABORT: the eight tubes are NOT equivalent (spread {spread:.2f} ps); "
                     "sigma_e must then be computed per tube, not once.")
sgm = float(sg.mean())
se = float(np.hypot(sgm, SIGMA_T_PAIR))
print(f"\n  Sim-1 sigma_geo   = {1000*sgm:6.1f} ps   (mean over the eight)")
print(f"  Sim-1 sigma_e     = {1000*se:6.1f} ps   (quadrature with "
      f"sigma_t/sqrt(2) = {1000*SIGMA_T_PAIR:.1f} ps)")
NEMA_GEO, NEMA_E, MERC_E, SHIPPED = 351.4, 357.2, 312.4, 312.4
print(f"\n  for comparison   NEMA sigma_geo {NEMA_GEO:.1f} ps, sigma_e {NEMA_E:.1f} ps")
print(f"                   Mercolli/shipped sigma_e {SHIPPED:.1f} ps")
ok = 1000 * sgm < NEMA_GEO
print(f"\n  REGISTERED: Sim-1 sigma_geo < NEMA's {NEMA_GEO:.1f} ps  ->  "
      f"{'CONFIRMED' if ok else 'REFUTED -- the geometry argument is wrong'}")
print(f"  Using the shipped 312.4 ps for Sim-1 would be wrong by "
      f"{100*(se-SHIPPED/1000)/(SHIPPED/1000):+.1f} %")
out.update(sigma_geo_ps=1000 * sgm, sigma_e_ps=1000 * se,
           sigma_e_ns=se, sigma_t_pair_ps=1000 * SIGMA_T_PAIR,
           symmetry_spread_ps=spread, nema_sigma_e_ps=NEMA_E,
           shipped_sigma_e_ps=SHIPPED, registered_smaller_than_nema=bool(ok))
json.dump(out, open(os.path.join(os.environ.get("PLI_RES") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                 "2026-08-26_A2_cond2/sigma_geo_sim1.json"), "w"), indent=1)
