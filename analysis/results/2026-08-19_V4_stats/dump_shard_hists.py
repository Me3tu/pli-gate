#!/usr/bin/env python3
"""V4(d) support — per-shard delay histograms for BOTH matching rules.

The two rules act on the SAME events, so tau_closest and tau_all are correlated
and sigma(Delta tau) computed as hypot(sigma_c, sigma_a) is conservative. To
measure the correlation we need paired resampling, which needs the per-shard
contributions rather than the pooled totals the production extractor keeps.

This emits, per VOI and per arm, the signal and delayed-window histograms
SHARD BY SHARD, on exactly the production TAU_BINS. Bootstrapping shards then
resamples both arms coherently, because a shard drawn for one arm is the same
shard for the other.

Physics untouched: every function is imported from the production extractor.
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
import h5py, numpy as np

EX = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-12_steinberger_allprompts"
sys.path.insert(0, str(EX))
import steinberger_extract_ap as X  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--singles", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--phantom", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--sigma-t-ps", type=float, default=90.9)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    cfg = json.load(open(a.config))[a.phantom]
    vois, outline = cfg["vois"], cfg["phantom_outline"]
    nb = len(X.TAU_BINS) - 1
    keys = [-1] + list(range(len(vois)))

    f = h5py.File(a.singles, "r"); g = f["singles"]
    sh = g["shard_id"][:]; bounds = X.find_shard_boundaries(sh); del sh
    ns = len(bounds)
    print(f"[{a.label}] {ns} shards, {len(vois)} VOIs", flush=True)

    H = {m: {w: np.zeros((ns, len(keys), nb)) for w in ("sig", "rnd")}
         for m in ("closest", "all")}
    for m in ("closest", "all"):
        rng = np.random.default_rng(int(a.sigma_t_ps * 1000) or 1)   # same seed per arm, as production
        for si, (sid, s, e) in enumerate(bounds):
            out = X.clinical_per_shard(
                g["time_ns"][s:e], g["in_511"][s:e], g["in_prompt"][s:e],
                g["gx"][s:e], g["gy"][s:e], g["gz"][s:e],
                a.sigma_t_ps / 1000.0, rng, vois, outline, matching=m)
            if out is None:
                continue
            by_sig, by_rand = out
            for ki, k in enumerate(keys):
                if by_sig[k]:
                    H[m]["sig"][si, ki] = np.histogram(np.concatenate(by_sig[k]), bins=X.TAU_BINS)[0]
                if by_rand[k]:
                    H[m]["rnd"][si, ki] = np.histogram(np.concatenate(by_rand[k]), bins=X.TAU_BINS)[0]
            if si % 32 == 0:
                print(f"  {m}: shard {si}/{ns}", flush=True)
    f.close()

    np.savez_compressed(a.out,
        tau_bins=X.TAU_BINS,
        voi_names=np.array(["bg"] + [v["name"] for v in vois]),
        tau_cfg=np.array([np.nan] + [float(v.get("tau_cfg_ns", np.nan)) for v in vois]),
        closest_sig=H["closest"]["sig"], closest_rnd=H["closest"]["rnd"],
        all_sig=H["all"]["sig"], all_rnd=H["all"]["rnd"])
    print(f"[{a.label}] wrote {a.out}")


if __name__ == "__main__":
    main()
