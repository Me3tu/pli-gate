#!/usr/bin/env python3
"""Operating point of the windowed-moment readout on the corrected-yield eight-tube acquisition (C0, zero delay,
production builder adj2, all-prompts matching): per tube, in the fit window (1, 5) ns of path-corrected delay,
the delayed-arm share of the prompt-arm counts (f_R, the definition of Note S11) and the subtracted signal
count N_s. Source: 2026-10-09_pair_window/c0/hists_pair_window.npz. Writes moment_window_point.json."""
import json, os
from pathlib import Path
import numpy as np
D = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-10-09_pair_window/c0"
z = np.load(D / "hists_pair_window.npz"); bins = z["bins"]
bc = 0.5 * (bins[:-1] + bins[1:]); fw = (bc > 1.0) & (bc < 5.0)
names = [r["voi"] for r in json.load(open(D / "pair_window.json"))["rows"]]
out = dict(source=str(D / "hists_pair_window.npz"), arm="adj2", window_ns=[1.0, 5.0], tubes={})
for rule in ("all", "closest"):
    hs = z[f"adj2_{rule}_sig"].sum(0)[:, fw].sum(1); hr = z[f"adj2_{rule}_rnd"].sum(0)[:, fw].sum(1)
    for v, n in enumerate(names, start=1):
        out["tubes"].setdefault(n, {})[rule] = dict(prompt=float(hs[v]), delayed=float(hr[v]), fR=float(hr[v] / hs[v]), Ns=float(hs[v] - hr[v]))
fr = [t["all"]["fR"] for t in out["tubes"].values()]; ns = [t["all"]["Ns"] for t in out["tubes"].values()]
out.update(fR_range=[min(fr), max(fr)], fR_median=float(np.median(fr)), Ns_range=[min(ns), max(ns)], Ns_median=float(np.median(ns)))
json.dump(out, open(Path(__file__).resolve().parent / "moment_window_point.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "tubes"}, indent=1))
for n, t in out["tubes"].items(): print(n, {r: (round(x["fR"], 3), int(x["Ns"])) for r, x in t.items()})
