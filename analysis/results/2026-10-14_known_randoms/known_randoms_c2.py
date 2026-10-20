#!/usr/bin/env python3
"""Review round 14 (2026-10-14): the known-content test repeated on the realistic-timing arm C2 (effective 388 ps
direct component, identical to C0 except for the direct delay). Copy of known_randoms_c0.py with SRC and outputs
changed. The validation against the C0 ladder is replaced by a count comparison of the k = 1 rung with the C2 production
histograms (2026-08-26_B2_paired/hists_sim1_388.npz, independent smear realization, so not bin for bin).

Original docstring (round 13): Review round 13 (2026-10-14), referee M4: a direct test of the delayed-window randoms correction against
known accidental content, at twelve accidental rates, on the corrected-yield zero-delay eight-tube stream C0.

Construction. The rate replay of 2026-10-08_c0_ladder (ladder_joint_c0.py) is repeated with the published
clinical extractor logic (sim1_extract.clinical_per_shard, re-implemented here line by line so that each triple
keeps the indices of its three singles, and the signal and randoms histograms are checked bin for bin against
ladder_joint_hists.npz). Compressing the decay epochs (t_new = t0/k + (t - t0) per eventID group) injects
accidental prompts and accidental pairs at a known rate while every decay keeps its own internal timing, so the
content of the signal window is known for every triple from the event identity of its three singles, which is
used ONLY to label, never to select:
  G        genuine: pair and prompt from one decay
  AP       accidental prompt on a genuine pair (pair from one decay, prompt from another)
  PA_corr  accidental pair whose prompt belongs to one of its two photons (time-correlated with the pair)
  PA_unc   accidental pair with an unrelated prompt
The delayed window (prompt search shifted by 50 ns) is meant to estimate the uncorrelated part, AP + PA_unc.
Histograms per shard, rung, rule, region: signal, randoms arm, and the four labeled classes of the signal arm.
Smear draws per shard rng [90900, shard], identical for every rung and rule (as the ladder).
Output known_randoms_hists.npz: hist[shard, k, rule, voi+1, cls, bin] with cls = sig, rnd, G, AP, PA_corr, PA_unc.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import h5py

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import sim1_extract as SX  # noqa: E402

SRC = str(Path(os.environ.get("PLI_RAW", "<missing>")) / "c2_sim1_388/singles_c2_sim1_388.h5")
CFG = json.loads((RES / "2026-08-17_derivation_test/sim1_voi_config.json").read_text())["sim1_tubes"]
VOIS, OUTLINE = CFG["vois"], CFG["phantom_outline"]
KS = [0.15, 0.25, 0.35, 0.50, 0.70, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00]
RULES = ["closest", "all"]
CLS = ["sig", "rnd", "G", "AP", "PA_corr", "PA_unc"]
SIGMA_T_NS = 0.0909
NB = len(SX.TAU_BINS) - 1
OUT = RES / "2026-10-14_known_randoms"


def rescale_times(t, eid, k):                     # copied from ladder_joint_c0.py
    order = np.argsort(eid, kind="stable")
    es, ts = eid[order], t[order]
    newgrp = np.empty(len(es), dtype=bool); newgrp[0] = True; newgrp[1:] = es[1:] != es[:-1]
    gid = np.cumsum(newgrp) - 1
    base = np.minimum.reduceat(ts, np.flatnonzero(newgrp))[gid]
    out = np.empty_like(t)
    out[order] = base / k + (ts - base)
    return out


def clinical_labeled(t_true, eid, in_511, in_prompt, gx, gy, gz, sigma_t_ns, rng, matching):
    """sim1_extract.clinical_per_shard, line by line, returning per-region histograms of the signal arm, the
    randoms arm, and the labeled classes of the signal arm."""
    H = np.zeros((len(VOIS) + 1, len(CLS), NB), dtype=np.int32)
    n = len(t_true)
    if n == 0: return H
    t_smear = t_true + rng.normal(0.0, sigma_t_ns, size=n)
    in_511_b = in_511.astype(bool); in_prompt_b = in_prompt.astype(bool)
    idx_511 = np.flatnonzero(in_511_b); idx_prompt = np.flatnonzero(in_prompt_b)
    if len(idx_511) < 2 or len(idx_prompt) < 1: return H
    t_511 = t_smear[idx_511]; t_prompt = t_smear[idx_prompt]
    o5 = np.argsort(t_511, kind="stable"); op = np.argsort(t_prompt, kind="stable")
    idx_511 = idx_511[o5]; t_511 = t_511[o5]
    idx_prompt = idx_prompt[op]; t_prompt = t_prompt[op]
    dt = np.diff(t_511); pair_ok = dt <= SX.W_PAIR_NS
    if not pair_ok.any(): return H
    pa = idx_511[:-1][pair_ok]; pb = idx_511[1:][pair_ok]
    p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
    p_ann = SX.tof_p_ann(p1, p2, t_smear[pa], t_smear[pb])
    perp = SX.lor_perp(p_ann, p1, p2)
    geom = (perp <= OUTLINE.get("lor_tol_mm", 10.0)) & SX.in_phantom_outline(p_ann[:, 0], p_ann[:, 1], p_ann[:, 2], OUTLINE)
    if not geom.any(): return H
    pa = pa[geom]; pb = pb[geom]; p_ann = p_ann[geom]
    t_pair_mean = 0.5 * (t_smear[pa] + t_smear[pb])
    if matching == "all":
        sel_s, pick_s = SX.all_prompt_indices(t_pair_mean, t_prompt, SX.W_PROMPT_NS)
        sel_r, pick_r = SX.all_prompt_indices(t_pair_mean + SX.DELAY_NS, t_prompt, SX.W_PROMPT_NS)
    else:
        has_s, pick_s = SX.closest_prompt_indices(t_pair_mean, t_prompt, SX.W_PROMPT_NS)
        has_r, pick_r = SX.closest_prompt_indices(t_pair_mean + SX.DELAY_NS, t_prompt, SX.W_PROMPT_NS)
        sel_s = np.flatnonzero(has_s); sel_r = np.flatnonzero(has_r)

    def build(sel, pick):
        if len(sel) == 0: return None, None, None
        pa_h = pa[sel]; pb_h = pb[sel]; p_ann_h = p_ann[sel]
        ip = idx_prompt[pick]
        p1_h = np.stack([gx[pa_h], gy[pa_h], gz[pa_h]], -1); p2_h = np.stack([gx[pb_h], gy[pb_h], gz[pb_h]], -1)
        pp_h = np.stack([gx[ip], gy[ip], gz[ip]], -1)
        tau = SX.per_event_tau(p1_h, p2_h, pp_h, t_smear[pa_h], t_smear[pb_h], t_smear[ip], p_ann_h)
        ea, eb, ep = eid[pa_h], eid[pb_h], eid[ip]
        pair_g = ea == eb
        lab = np.where(pair_g & (ep == ea), 2, np.where(pair_g, 3, np.where((ep == ea) | (ep == eb), 4, 5)))
        return tau, p_ann_h, lab

    tau_s, p_ann_s, lab_s = build(sel_s, pick_s)
    tau_r, p_ann_r, _ = build(sel_r, pick_r)
    if tau_r is not None: tau_r = tau_r + SX.DELAY_NS
    if tau_s is not None:
        ts32 = tau_s.astype(np.float32)
        masks = SX.voi_masks(p_ann_s, VOIS) + [SX.voi_bg_mask(p_ann_s, VOIS)]
        for i, m in enumerate(masks):
            vi = i + 1 if i < len(VOIS) else 0
            if not m.any(): continue
            H[vi, 0] = np.histogram(ts32[m], bins=SX.TAU_BINS)[0]
            for c in (2, 3, 4, 5):
                mc = m & (lab_s == c)
                if mc.any(): H[vi, c] = np.histogram(ts32[mc], bins=SX.TAU_BINS)[0]
    if tau_r is not None:
        tr32 = tau_r.astype(np.float32)
        masks = SX.voi_masks(p_ann_r, VOIS) + [SX.voi_bg_mask(p_ann_r, VOIS)]
        for i, m in enumerate(masks):
            vi = i + 1 if i < len(VOIS) else 0
            if m.any(): H[vi, 1] = np.histogram(tr32[m], bins=SX.TAU_BINS)[0]
    return H


def one_shard(args):
    sh_id, a, b, src = args
    with h5py.File(src, "r") as f:
        g = f["singles"]
        t = g["time_ns"][a:b].astype(np.float64); eid = g["eventID"][a:b]
        in5 = g["in_511"][a:b]; inp = g["in_prompt"][a:b]
        gx = g["gx"][a:b]; gy = g["gy"][a:b]; gz = g["gz"][a:b]
    H = np.zeros((len(KS), len(RULES), len(VOIS) + 1, len(CLS), NB), dtype=np.int32)
    for ik, k in enumerate(KS):
        tk = t if k == 1.0 else rescale_times(t, eid, k)
        for ir, rule in enumerate(RULES):
            rng = np.random.default_rng([90900, int(sh_id)])
            H[ik, ir] = clinical_labeled(tk, eid, in5, inp, gx, gy, gz, SIGMA_T_NS, rng, rule)
    return sh_id, H


def main():
    if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
    OUT.mkdir(parents=True, exist_ok=True)
    with h5py.File(SRC, "r") as f:
        sh = f["singles"]["shard_id"][:]
    bounds = SX.find_shard_boundaries(sh)
    print(f"{len(bounds)} shards, {len(KS)} rungs, {len(RULES)} rules", flush=True)
    t0 = time.time(); res = {}
    with Pool(24) as pool:
        for i, (sid, H) in enumerate(pool.imap_unordered(one_shard, [(s, a, b, SRC) for s, a, b in bounds])):
            res[sid] = H
            if i % 8 == 0: print(f"  {i+1}/{len(bounds)} shards  {time.time()-t0:.0f} s", flush=True)
    ids = sorted(res)
    np.savez_compressed(OUT / "known_randoms_c2_hists.npz", hist=np.stack([res[s] for s in ids]), shard_ids=np.array(ids),
                        ks=np.array(KS), rules=np.array(RULES), classes=np.array(CLS), bins=SX.TAU_BINS,
                        voi_names=np.array(["bg"] + [v["name"] for v in VOIS]),
                        tau_cfg=np.array([np.nan] + [v["tau_cfg_ns"] for v in VOIS]))
    # validation: labels partition the signal arm, and k = 1 counts against the C2 production histograms
    HK = np.stack([res[s] for s in ids]); lab_sum = bool(np.array_equal(HK[..., 2:, :].sum(-2), HK[..., 0, :]))
    P = np.load(RES / "2026-08-26_B2_paired/hists_sim1_388.npz", allow_pickle=True); ik = KS.index(1.00); cmp = {}
    for r, rule in enumerate(RULES):
        for c, cls in ((0, "sig"), (1, "rnd")):
            a = float(HK[:, ik, r, :, c, :].sum()); b = float(P[f"{rule}_{cls}"].sum()); cmp[f"{rule}_{cls}"] = [a, b, (a - b) / b]
    print(f"validation: labels partition the signal arm: {lab_sum}; k=1 counts vs production: {cmp}")
    json.dump(dict(labels_partition_signal=lab_sum, k1_counts_vs_production=cmp, wall_s=time.time() - t0),
              open(OUT / "known_randoms_c2_validation.json", "w"), indent=1)
    print(f"done {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
