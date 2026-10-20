#!/usr/bin/env python3
"""Review round 16 (2026-10-19): prompt-selection rules in use, the all-pairs builder, and the pair-breaking
mechanism, with known-content labels, on the corrected-yield zero-delay eight-tube stream C0.

Base: 2026-10-14_known_randoms/known_randoms_c0.py (clinical extractor of sim1_extract, re-implemented so that each
triple keeps the indices of its singles; event identity labels triples and never selects). Rate replay rungs
k = 0.5, 1.0, 2.0 (k = 1.0 is the production rate, mu = 0.029). Smear draws per shard rng [90900, shard], identical
for every rung and arm, as in known_randoms_c0.py.

Arms (signal window [T-5, T+5] ns around the pair time T, randoms arm the same rule at T + 50 ns):
  closest      production closest-prompt rule (sim1_extract.closest_prompt_indices; validation)
  all          all-prompts matching (validation against known_randoms_hists.npz at k = 1.0)
  closest_x    exact nearest prompt among all candidates (ties to the earlier), to check the production rule,
               which compares only the earliest and the latest candidate
  first        first-arriving prompt in the window (list-mode 'first qualifying' rule)
  highest      highest-energy prompt in the window
  veto1        the pair is kept only if exactly one prompt lies in the window, in both windows (single-candidate
               veto applied identically in the delayed window)
  ap_all       all-pairs builder (every pair of annihilation singles within 2 ns, Note S18) with all-prompts matching
  ap_closest   all-pairs builder with the production closest-prompt rule
Labels: G genuine, AP accidental prompt on a genuine pair, PA_corr accidental pair carrying its own prompt, PA_unc.
Also, at k = 1.0, a truth diagnostic of pair breaking: for every decay whose two annihilation singles fall within
2 ns, whether the adjacent-singles builder forms that pair, split by the presence of another decay's prompt-window
single in [T-5, T+5] (signal) and in [T+45, T+55] (delayed).
Output: rules16_hists.npz, rules16_diag.json."""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import h5py

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
sys.path.insert(0, str(RES / "2026-08-17_derivation_test"))
import sim1_extract as SX  # noqa: E402

SRC = str(Path(os.environ.get("PLI_RAW", "<missing>")) / "feedback3_c0/c0_sim1_0ps/singles_c0_sim1_0ps.h5")
CFG = json.loads((RES / "2026-08-17_derivation_test/sim1_voi_config.json").read_text())["sim1_tubes"]
VOIS, OUTLINE = CFG["vois"], CFG["phantom_outline"]
KS = [0.5, 1.0, 2.0]
ARMS = ["closest", "all", "closest_x", "first", "highest", "veto1", "ap_all", "ap_closest"]
CLS = ["sig", "rnd", "G", "AP", "PA_corr", "PA_unc"]
SIGMA_T_NS = 0.0909
NB = len(SX.TAU_BINS) - 1
OUT = RES / "2026-10-19_review16"


def rescale_times(t, eid, k):                     # as known_randoms_c0.py
    order = np.argsort(eid, kind="stable")
    es, ts = eid[order], t[order]
    newgrp = np.empty(len(es), dtype=bool); newgrp[0] = True; newgrp[1:] = es[1:] != es[:-1]
    gid = np.cumsum(newgrp) - 1
    base = np.minimum.reduceat(ts, np.flatnonzero(newgrp))[gid]
    out = np.empty_like(t)
    out[order] = base / k + (ts - base)
    return out


def adjacent_pairs(n):
    return np.arange(n - 1), np.arange(1, n)


def allpairs(t_511, w):
    A, B = [], []
    k = 1
    while k < len(t_511):
        ok = (t_511[k:] - t_511[:-k]) <= w
        if not ok.any(): break
        ia = np.flatnonzero(ok); A.append(ia); B.append(ia + k); k += 1
    return (np.concatenate(A), np.concatenate(B)) if A else (np.zeros(0, int), np.zeros(0, int))


def select(rule, T, t_prompt, e_prompt):
    """Return (pair_index, prompt_index) of the triples formed for pair times T."""
    W = SX.W_PROMPT_NS
    if rule == "closest":
        has, pick = SX.closest_prompt_indices(T, t_prompt, W)
        return np.flatnonzero(has), pick
    if rule == "all":
        return SX.all_prompt_indices(T, t_prompt, W)
    lo = np.searchsorted(t_prompt, T - W, side="left"); hi = np.searchsorted(t_prompt, T + W, side="right")
    n = hi - lo; has = n > 0
    if rule == "first":
        return np.flatnonzero(has), lo[has]
    if rule == "veto1":
        one = n == 1
        return np.flatnonzero(one), lo[one]
    pi, pk = SX.all_prompt_indices(T, t_prompt, W)
    if len(pi) == 0: return pi, pk
    if rule == "closest_x":
        key = np.abs(t_prompt[pk] - T[pi])
        o = np.lexsort((t_prompt[pk], key, pi))
    elif rule == "highest":
        o = np.lexsort((t_prompt[pk], -e_prompt[pk], pi))
    else:
        raise ValueError(rule)
    pi, pk = pi[o], pk[o]
    first = np.ones(len(pi), bool); first[1:] = pi[1:] != pi[:-1]
    return pi[first], pk[first]


def arm_hists(arm, t_smear, eid, idx_511, t_511, idx_prompt, t_prompt, e_prompt, gx, gy, gz):
    H = np.zeros((len(VOIS) + 1, len(CLS), NB), dtype=np.int32)
    builder = "allpairs" if arm.startswith("ap_") else "adjacent"
    rule = arm[3:] if arm.startswith("ap_") else arm
    if builder == "adjacent":
        dt = np.diff(t_511); ok = dt <= SX.W_PAIR_NS
        ia = np.flatnonzero(ok); ib = ia + 1
    else:
        ia, ib = allpairs(t_511, SX.W_PAIR_NS)
    if len(ia) == 0: return H
    pa = idx_511[ia]; pb = idx_511[ib]
    p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1); p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
    p_ann = SX.tof_p_ann(p1, p2, t_smear[pa], t_smear[pb])
    perp = SX.lor_perp(p_ann, p1, p2)
    geom = (perp <= OUTLINE.get("lor_tol_mm", 10.0)) & SX.in_phantom_outline(p_ann[:, 0], p_ann[:, 1], p_ann[:, 2], OUTLINE)
    if not geom.any(): return H
    pa = pa[geom]; pb = pb[geom]; p_ann = p_ann[geom]
    T = 0.5 * (t_smear[pa] + t_smear[pb])
    sel_s, pick_s = select(rule, T, t_prompt, e_prompt)
    sel_r, pick_r = select(rule, T + SX.DELAY_NS, t_prompt, e_prompt)

    def build(sel, pick):
        if len(sel) == 0: return None, None, None
        pa_h = pa[sel]; pb_h = pb[sel]; p_ann_h = p_ann[sel]; ip = idx_prompt[pick]
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


def pair_break_diag(t_smear, eid, idx_511, t_511, idx_prompt, t_prompt):
    """Truth diagnostic at one rung: for decays with exactly two annihilation-window singles within 2 ns, is the pair
    formed by the adjacent builder, by presence of another decay's prompt-window single near T and near T + 50 ns."""
    pos = np.empty(len(idx_511), np.int64); pos[:] = np.arange(len(idx_511))
    e5 = eid[idx_511]
    o = np.lexsort((t_511, e5)); es = e5[o]
    start = np.ones(len(es), bool); start[1:] = es[1:] != es[:-1]
    gstart = np.flatnonzero(start); gsize = np.diff(np.append(gstart, len(es)))
    two = gstart[gsize == 2]
    a, b = o[two], o[two + 1]                      # positions in the time-sorted 511 list
    dt = np.abs(t_511[b] - t_511[a]); ok = dt <= SX.W_PAIR_NS
    a, b = a[ok], b[ok]
    lo_, hi_ = np.minimum(a, b), np.maximum(a, b)
    formed = (hi_ - lo_) == 1
    T = 0.5 * (t_511[a] + t_511[b]); ev = e5[a]
    ep = eid[idx_prompt]

    def n_other(Tc):
        l = np.searchsorted(t_prompt, Tc - SX.W_PROMPT_NS, "left"); h = np.searchsorted(t_prompt, Tc + SX.W_PROMPT_NS, "right")
        cnt = h - l
        # subtract own prompt(s) of the same decay inside the window
        pi = np.repeat(np.arange(len(Tc)), cnt)
        pk = np.concatenate([np.arange(x, y) for x, y in zip(l, h)]) if cnt.sum() else np.zeros(0, int)
        own = np.bincount(pi[ep[pk] == ev[pi]], minlength=len(Tc)) if len(pi) else np.zeros(len(Tc), int)
        return cnt - own

    ns, nd = n_other(T), n_other(T + SX.DELAY_NS)
    out = dict(n_pairs=int(len(T)), n_formed=int(formed.sum()))
    for name, m in (("sig_acc", ns > 0), ("del_acc", nd > 0), ("sig_none", ns == 0), ("del_none", nd == 0)):
        out[name] = [int(m.sum()), int((formed & m).sum())]
    return out


def multiplicity_counts(t_511, t_prompt):
    dt = np.diff(t_511); ia = np.flatnonzero(dt <= SX.W_PAIR_NS)
    T = 0.5 * (t_511[ia] + t_511[ia + 1])
    lo = np.searchsorted(t_prompt, T - SX.W_PROMPT_NS, "left"); hi = np.searchsorted(t_prompt, T + SX.W_PROMPT_NS, "right")
    n = hi - lo
    return [int((n == j).sum()) for j in range(4)] + [int((n >= 4).sum())]


def one_shard(args):
    sh_id, a, b, src = args
    with h5py.File(src, "r") as f:
        g = f["singles"]
        t = g["time_ns"][a:b].astype(np.float64); eid = g["eventID"][a:b]
        in5 = g["in_511"][a:b]; inp = g["in_prompt"][a:b]; en = g["energy_keV"][a:b]
        gx = g["gx"][a:b]; gy = g["gy"][a:b]; gz = g["gz"][a:b]
    H = np.zeros((len(KS), len(ARMS), len(VOIS) + 1, len(CLS), NB), dtype=np.int32)
    diag, mult = None, None
    for ik, k in enumerate(KS):
        tk = t if k == 1.0 else rescale_times(t, eid, k)
        rng = np.random.default_rng([90900, int(sh_id)])
        t_smear = tk + rng.normal(0.0, SIGMA_T_NS, size=len(tk))
        idx_511 = np.flatnonzero(in5.astype(bool)); idx_prompt = np.flatnonzero(inp.astype(bool))
        if len(idx_511) < 2 or len(idx_prompt) < 1: continue
        t_511 = t_smear[idx_511]; t_prompt = t_smear[idx_prompt]
        o5 = np.argsort(t_511, kind="stable"); op = np.argsort(t_prompt, kind="stable")
        idx_511 = idx_511[o5]; t_511 = t_511[o5]; idx_prompt = idx_prompt[op]; t_prompt = t_prompt[op]
        e_prompt = en[idx_prompt].astype(np.float64)
        for ia_, arm in enumerate(ARMS):
            H[ik, ia_] = arm_hists(arm, t_smear, eid, idx_511, t_511, idx_prompt, t_prompt, e_prompt, gx, gy, gz)
        if k == 1.0:
            diag = pair_break_diag(t_smear, eid, idx_511, t_511, idx_prompt, t_prompt)
            mult = multiplicity_counts(t_511, t_prompt)
    return sh_id, H, diag, mult


def main():
    if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
    OUT.mkdir(parents=True, exist_ok=True)
    with h5py.File(SRC, "r") as f:
        sh = f["singles"]["shard_id"][:]
    bounds = SX.find_shard_boundaries(sh)
    print(f"{len(bounds)} shards, rungs {KS}, arms {ARMS}", flush=True)
    t0 = time.time(); res, diags, mults = {}, [], []
    with Pool(20) as pool:
        for i, (sid, H, d, m) in enumerate(pool.imap_unordered(one_shard, [(s, a, b, SRC) for s, a, b in bounds])):
            res[sid] = H
            if d: diags.append(d)
            if m: mults.append(m)
            if i % 8 == 0: print(f"  {i+1}/{len(bounds)} shards  {time.time()-t0:.0f} s", flush=True)
    ids = sorted(res); HK = np.stack([res[s] for s in ids])
    np.savez_compressed(OUT / "rules16_hists.npz", hist=HK, shard_ids=np.array(ids), ks=np.array(KS),
                        arms=np.array(ARMS), classes=np.array(CLS), bins=SX.TAU_BINS,
                        voi_names=np.array(["bg"] + [v["name"] for v in VOIS]),
                        tau_cfg=np.array([np.nan] + [v["tau_cfg_ns"] for v in VOIS]))
    D = {k: (sum(d[k] for d in diags) if isinstance(diags[0][k], int) else [sum(d[k][j] for d in diags) for j in range(2)])
         for k in diags[0]}
    M = [sum(m[j] for m in mults) for j in range(5)]
    # validation against known_randoms_hists.npz at k = 1.0 (index 5 there)
    KR = np.load(RES / "2026-10-14_known_randoms/known_randoms_hists.npz")
    HR = KR["hist"][:, list(KR["ks"]).index(1.0)]
    same = {arm: bool(np.array_equal(HK[:, KS.index(1.0), ARMS.index(arm)], HR[:, list(KR["rules"]).index(arm)]))
            for arm in ("closest", "all")}
    lab_sum = bool(np.array_equal(HK[..., 2:, :].sum(-2), HK[..., 0, :]))
    json.dump(dict(doc=__doc__, identical_to_known_randoms=same, labels_partition_signal=lab_sum,
                   pair_break=D, prompt_multiplicity_k1=dict(counts_0_1_2_3_ge4=M), wall_s=time.time() - t0),
              open(OUT / "rules16_diag.json", "w"), indent=1)
    print("validation", same, "labels partition", lab_sum, flush=True)
    print("pair_break", D, flush=True); print("multiplicity", M, flush=True)
    print(f"done {time.time()-t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
