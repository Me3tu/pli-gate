#!/usr/bin/env python3
"""mu = mean number of prompts inside a +/-5 ns window, per configuration.

Computed PER SHARD and then averaged. Computing across shard boundaries would
be wrong: shards are independent time streams, so a gap spanning a boundary is
not a real gap.

The prompt definition is each file's OWN stored in_prompt flag, so the isotope
window is whatever that dataset used (Sc-44 [1100,1200], I-124 [568,639]).
"""
from __future__ import annotations
import json, os, sys
from pathlib import Path
import numpy as np
import h5py

W_PROMPT_NS = 5.0
RAW = Path(os.environ.get("PLI_RAW", "<missing>"))
FILES = [
    ("Sim-1 frozen 0 ps", str(RAW / "2026-06-10_brain_tumor_1B_truetime/sim1_singles_truetime.h5")),
    ("Sim-1 C2 low-rate 388 ps", str(RAW / "c2_sim1_388_LOWRATE_null/singles_c2_sim1_388.h5")),
    ("NEMA T6 0 ps", str(RAW / "2026-08-16_T6_nema/singles_nema_all.h5")),
    ("NEMA C1 388 ps", str(RAW / "c1_nema388/singles_c1_nema388.h5")),
    ("Steinberger disks", str(RAW / "pli-quadra-archive-2026-08-15-singles/singles_steinberger_disks.h5")),
    ("Steinberger water", str(RAW / "pli-quadra-archive-2026-08-15-singles/singles_steinberger_water.h5")),
]
if not any(Path(p).exists() for _, p in FILES): raise SystemExit(f"{RAW}: raw data not included in the release, see analysis/README.md")
MAX_SHARDS = 6          # enough for a stable mean. These files are 5-10 GB


def mu_for(path: str):
    p = Path(path)
    if not p.exists():
        return None
    with h5py.File(p, "r") as f:
        g = f["singles"] if "singles" in f else f
        n = g["time_ns"].shape[0]
        take = min(n, 40_000_000)
        sh = g["shard_id"][:take]
        t = g["time_ns"][:take]
        ip = g["in_prompt"][:take].astype(bool)
    mus, spans = [], []
    for s in np.unique(sh)[:MAX_SHARDS]:
        m = sh == s
        tp = np.sort(t[m & ip])
        if len(tp) < 1000:
            continue
        gap = float(np.mean(np.diff(tp)))
        mus.append(2.0 * W_PROMPT_NS / gap)
        spans.append(float(t[m].max() - t[m].min()))
    if not mus:
        return None
    return dict(mu=float(np.mean(mus)), mu_sd=float(np.std(mus, ddof=1)) if len(mus) > 1 else 0.0,
                n_shards=len(mus), span_ns=float(np.mean(spans)))


out = {}
print(f"{'configuration':28s} {'mu':>9s} {'+/-':>8s} {'shards':>7s} {'span (s)':>10s}")
for lab, path in FILES:
    r = mu_for(path)
    if r is None:
        print(f"{lab:28s} {'MISSING':>9s}"); continue
    out[lab] = r
    print(f"{lab:28s} {r['mu']:9.4f} {r['mu_sd']:8.4f} {r['n_shards']:7d} {r['span_ns']/1e9:10.4f}")

(Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) /
     "2026-08-26_rate_scan/mu_by_config.json").write_text(json.dumps(out, indent=1))
print("\n  wrote mu_by_config.json")
