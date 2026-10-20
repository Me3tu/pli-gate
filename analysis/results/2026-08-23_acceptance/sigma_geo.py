#!/usr/bin/env python3
"""sigma(Delta_geo) a priori, from the scanner geometry alone.

Delta_geo = [d_p - (d_1 + d_2)/2] / c, the path-length term that appears in the
box variable. It is fixed by the Quadra geometry and the phantom extent -- no
data and no fit enter. Detector: R = 392.83 to 414.26 mm, height 1060 mm, from
gate/geometry/quadra_detector.mac lines 40-46.
"""
from __future__ import annotations
import numpy as np

C = 299.792458                    # mm/ns
R_IN, R_OUT, H = 392.8322, 414.2571, 1060.0
R_DET = 0.5 * (R_IN + R_OUT)      # interaction radius, mid-crystal
HZ = 0.5 * H


def hit_distance(p, u):
    """Distance from p along unit direction u to the cylinder R_DET, |z| <= HZ."""
    a = u[:, 0] ** 2 + u[:, 1] ** 2
    b = 2 * (p[:, 0] * u[:, 0] + p[:, 1] * u[:, 1])
    c = p[:, 0] ** 2 + p[:, 1] ** 2 - R_DET ** 2
    disc = np.maximum(b * b - 4 * a * c, 0.0)
    t = (-b + np.sqrt(disc)) / (2 * a)
    z = p[:, 2] + t * u[:, 2]
    ok = np.abs(z) <= HZ           # inside the axial acceptance
    return t, ok


def run(n=2_000_000, phantom="mercolli", seed=20260823):
    rng = np.random.default_rng(seed)
    if phantom == "mercolli":      # four tubes, r = 6 mm, L = 45 mm, near axis
        r = 60.0 * np.sqrt(rng.random(n))      # tube ring radius spread
        th = rng.uniform(0, 2 * np.pi, n)
        p = np.c_[r * np.cos(th), r * np.sin(th), rng.uniform(-22.5, 22.5, n)]
    else:                          # water cylinder, R = 25 mm, L = 300 mm
        r = 25.0 * np.sqrt(rng.random(n)); th = rng.uniform(0, 2 * np.pi, n)
        p = np.c_[r * np.cos(th), r * np.sin(th), rng.uniform(-150, 150, n)]

    def iso(m):
        c = rng.uniform(-1, 1, m); s = np.sqrt(1 - c * c)
        ph = rng.uniform(0, 2 * np.pi, m)
        return np.c_[s * np.cos(ph), s * np.sin(ph), c]

    u = iso(n)                     # the 511 pair axis
    d1, ok1 = hit_distance(p, u)
    d2, ok2 = hit_distance(p, -u)
    v = iso(n)                     # the prompt, independent direction
    dp, okp = hit_distance(p, v)
    ok = ok1 & ok2 & okp
    delta = (dp[ok] - 0.5 * (d1[ok] + d2[ok])) / C
    return delta


for ph in ("mercolli", "water"):
    d = run(phantom=ph)
    print(f"{ph:9s}  n={len(d):,}  Delta_geo mean {1000*d.mean():+7.1f} ps   "
          f"sd {1000*d.std(ddof=1):6.1f} ps   "
          f"({d.std(ddof=1)*C:.0f} mm of path)")
SIG_T = 0.0909
print(f"\ntiming term  (eps1-eps2)/2 : sd = {1000*SIG_T/np.sqrt(2):.1f} ps")
d = run(phantom="water")
tot = np.sqrt(d.std(ddof=1) ** 2 + (SIG_T / np.sqrt(2)) ** 2)
print(f"quadrature total, water    : {1000*tot:.1f} ps")
print(f"fitted from the data (R1)  : 361.8 +/- 5.8 ps")
