#!/usr/bin/env python3
"""Event-identity (oracle) readout of the 20 cm water cylinder, as for the 20 mm cylinder (steinberger_extract_ap.oracle_run,
no smear, censored tail mean above 1 ns), per batch file and pooled with count weights. Writes w20_oracle.json."""
import json, os, sys
from pathlib import Path
import numpy as np
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[2])
sys.path.insert(0, str(RES / "2026-08-12_steinberger_allprompts")); sys.path.insert(0, str(Path(__file__).resolve().parent))
import steinberger_extract_ap as X  # noqa
import importlib.util
spec = importlib.util.spec_from_file_location("w20", str(Path(__file__).resolve().parent / "w20_extract.py"))
VOIS = [dict(name=n, shape="cylinder", cx_mm=0.0, cy_mm=0.0, cz_mm=0.0, r_mm=r, half_z_mm=hz, tau_cfg_ns=1.84)
        for n, r, hz in (("VOI_5x15", 25.0, 75.0), ("VOI_10x20", 50.0, 100.0), ("VOI_15x25", 75.0, 125.0))]
acc = {}
for fn in sys.argv[1:]:
    for r in X.oracle_run(fn, VOIS, sigma_t_ps=0.0):
        a = acc.setdefault(r["voi"], [0.0, 0])
        if r["n_signal"]: a[0] += r["tau_hat_ns"] * r["n_signal"]; a[1] += r["n_signal"]
out = {}
for k, (s, n) in acc.items():
    t = s / n; out[k] = dict(tau_hat_ns=t, sigma_ns=t / np.sqrt(n), n_signal=n)
    print(f"{k:10s} tau {t:.4f} +- {t/np.sqrt(n):.4f}  n {n}")
json.dump(out, open(Path(__file__).resolve().parent / "w20_oracle.json", "w"), indent=1)
