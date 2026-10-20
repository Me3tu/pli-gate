"""mu per rung on C0: expected number of prompt singles in a +-5 ns window = 10 ns / mean spacing,
measured after the same epoch compression, per shard; reported as mean over all shards and over
the first six shards (the convention of 2026-08-26_rate_scan/mu_ladder.json)."""
import json, os, sys
from pathlib import Path
import numpy as np, h5py
RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-10-08_c0_ladder"))
from ladder_joint_c0 import rescale_times, KS, SRC
sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import sim1_extract as SX
if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
with h5py.File(SRC, "r") as f:
    g = f["singles"]; sh = g["shard_id"][:]
    bounds = SX.find_shard_boundaries(sh)
    res = {f"{k:.2f}": [] for k in KS}
    for (sid, a, b) in bounds:
        t = g["time_ns"][a:b].astype(np.float64); eid = g["eventID"][a:b]; ip = g["in_prompt"][a:b].astype(bool)
        for k in KS:
            tk = t if k == 1.0 else rescale_times(t, eid, k)
            tp = np.sort(tk[ip])
            if len(tp) > 1000:
                res[f"{k:.2f}"].append(10.0 / float(np.mean(np.diff(tp))))
out = {kk: dict(k=float(kk), mu=float(np.mean(v)), mu_first6=float(np.mean(v[:6])), n_shards=len(v),
                sd_over_shards=float(np.std(v, ddof=1))) for kk, v in res.items()}
(RES / "2026-10-08_c0_ladder/mu_ladder_c0.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
