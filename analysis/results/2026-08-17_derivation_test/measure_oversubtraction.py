#!/usr/bin/env python3
"""N2 — measure the delayed-window over-subtraction directly, from truth tags.

No modelling. For every ROI we count, inside the fit window (tau_min < tau <
W_PROMPT):

    n_total      accepted triples in the signal arm
    n_sub        what the delayed-window estimator subtracts  (the randoms arm)
    n_true_acc   accepted signal-arm triples that are NOT genuine, by TRUTH
    n_true_sig   accepted signal-arm triples that ARE genuine, by TRUTH

A triple is genuine iff its two 511 singles and its prompt single all carry the
same eventID within the shard, i.e. all three came from one decay. eventID is
used ONLY to label, never to select -- the pipeline that produces the triples is
byte-for-byte the clinical one.

    over-subtraction  =  n_sub - n_true_acc
    phi               = (n_sub - n_true_acc) / n_true_sig

n_true_sig above tau_min = 1 ns is the o-Ps count: p-Ps (0.125 ns) and direct
annihilation (0.388 ns) are sub-nanosecond and contribute negligibly there.

Run under BOTH matching rules. Registered expectation (before running): under
all-prompts n_sub and n_true_acc should agree within counting error; under
closest-prompt n_sub should EXCEED n_true_acc.

Every physics function is imported from steinberger_extract_ap.py so the
estimator, windows and geometry are identical to the production runs.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

EXDIR = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-12_steinberger_allprompts"
sys.path.insert(0, str(EXDIR))
import steinberger_extract_ap as EX  # noqa: E402


def clinical_per_shard_truth(eventID, t_true, in_511, in_prompt, gx, gy, gz,
                             sigma_t_ns, rng, vois, phantom_outline, matching):
    """EX.clinical_per_shard with a truth flag carried alongside each triple."""
    n = len(t_true)
    if n == 0:
        return None
    t_smear = t_true + rng.normal(0.0, sigma_t_ns, size=n)
    in_511_b = in_511.astype(bool); in_prompt_b = in_prompt.astype(bool)
    idx_511 = np.flatnonzero(in_511_b); idx_prompt = np.flatnonzero(in_prompt_b)
    if len(idx_511) < 2 or len(idx_prompt) < 1:
        return None
    t_511 = t_smear[idx_511]; t_prompt = t_smear[idx_prompt]
    o5 = np.argsort(t_511); op = np.argsort(t_prompt)
    idx_511 = idx_511[o5]; t_511 = t_511[o5]
    idx_prompt = idx_prompt[op]; t_prompt = t_prompt[op]
    dt = np.diff(t_511); pair_ok = dt <= EX.W_PAIR_NS
    if not pair_ok.any():
        return None
    pa = idx_511[:-1][pair_ok]; pb = idx_511[1:][pair_ok]
    p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1)
    p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
    p_ann = EX.tof_p_ann(p1, p2, t_smear[pa], t_smear[pb])
    perp = EX.lor_perp(p_ann, p1, p2)
    geom = (perp <= phantom_outline.get("lor_tol_mm", 10.0)) & EX.in_phantom_outline(
        p_ann[:, 0], p_ann[:, 1], p_ann[:, 2], phantom_outline)
    if not geom.any():
        return None
    pa = pa[geom]; pb = pb[geom]; p_ann = p_ann[geom]
    t_pair_mean = 0.5 * (t_smear[pa] + t_smear[pb])

    if matching == "all":
        sel_s, pick_s = EX.all_prompt_indices(t_pair_mean, t_prompt, EX.W_PROMPT_NS)
        sel_r, pick_r = EX.all_prompt_indices(t_pair_mean + EX.DELAY_NS, t_prompt, EX.W_PROMPT_NS)
    else:
        has_s, pick_s = EX.closest_prompt_indices(t_pair_mean, t_prompt, EX.W_PROMPT_NS)
        has_r, pick_r = EX.closest_prompt_indices(t_pair_mean + EX.DELAY_NS, t_prompt, EX.W_PROMPT_NS)
        sel_s = np.flatnonzero(has_s); sel_r = np.flatnonzero(has_r)

    def build(sel, pick, want_truth):
        if len(sel) == 0:
            return None, None, None
        pa_h = pa[sel]; pb_h = pb[sel]; p_ann_h = p_ann[sel]
        ip = idx_prompt[pick]
        p1_h = np.stack([gx[pa_h], gy[pa_h], gz[pa_h]], -1)
        p2_h = np.stack([gx[pb_h], gy[pb_h], gz[pb_h]], -1)
        pp_h = np.stack([gx[ip], gy[ip], gz[ip]], -1)
        tau = EX.per_event_tau(p1_h, p2_h, pp_h,
                               t_smear[pa_h], t_smear[pb_h], t_smear[ip], p_ann_h)
        truth = None
        if want_truth:
            ea = eventID[pa_h]; eb = eventID[pb_h]; ep = eventID[ip]
            truth = (ea == eb) & (ea == ep)          # all three from one decay
        return tau, p_ann_h, truth

    tau_s, p_ann_s, tru_s = build(sel_s, pick_s, True)
    tau_r, p_ann_r, _ = build(sel_r, pick_r, False)
    if tau_r is not None:
        tau_r = tau_r + EX.DELAY_NS

    by_sig = {i: [] for i in range(-1, len(vois))}
    by_tru = {i: [] for i in range(-1, len(vois))}
    by_rand = {i: [] for i in range(-1, len(vois))}
    if tau_s is not None:
        ms = EX.voi_masks(p_ann_s, vois)
        for i, m in enumerate(ms):
            if m.any():
                by_sig[i].append(tau_s[m].astype(np.float32))
                by_tru[i].append(tru_s[m])
        b = EX.voi_bg_mask(p_ann_s, vois)
        if b.any():
            by_sig[-1].append(tau_s[b].astype(np.float32)); by_tru[-1].append(tru_s[b])
    if tau_r is not None:
        mr = EX.voi_masks(p_ann_r, vois)
        for i, m in enumerate(mr):
            if m.any():
                by_rand[i].append(tau_r[m].astype(np.float32))
        b = EX.voi_bg_mask(p_ann_r, vois)
        if b.any():
            by_rand[-1].append(tau_r[b].astype(np.float32))
    return by_sig, by_tru, by_rand


def counts_in_fit_window(tau_sig, truth_sig, tau_rand):
    """Counts inside exactly the window subtract_and_mle fits."""
    lo, hi = EX.TAU_MIN_NS, EX.W_PROMPT_NS
    s = (tau_sig > lo) & (tau_sig < hi)
    r = (tau_rand > lo) & (tau_rand < hi)
    n_total = int(s.sum())
    n_sub = int(r.sum())
    n_true_sig = int((s & truth_sig).sum())
    n_true_acc = int((s & ~truth_sig).sum())
    return n_total, n_sub, n_true_sig, n_true_acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--singles", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--phantom", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--matching", choices=["closest", "all"], required=True)
    ap.add_argument("--sigma-t-ps", type=float, default=90.9)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    cfg = json.load(open(args.config))[args.phantom]
    vois = cfg["vois"]; outline = cfg["phantom_outline"]
    sigma_t_ns = args.sigma_t_ps / 1000.0
    rng = np.random.default_rng(int(args.sigma_t_ps * 1000) or 1)

    f = h5py.File(args.singles, "r"); g = f["singles"]
    sh = g["shard_id"][:]; bounds = EX.find_shard_boundaries(sh); del sh
    print(f"[{args.label}/{args.matching}] {len(bounds)} shards, {len(vois)} VOIs", flush=True)
    acc_s = {i: [] for i in range(-1, len(vois))}
    acc_t = {i: [] for i in range(-1, len(vois))}
    acc_r = {i: [] for i in range(-1, len(vois))}
    for k, (sh_id, s, e) in enumerate(bounds):
        out = clinical_per_shard_truth(
            g["eventID"][s:e], g["time_ns"][s:e], g["in_511"][s:e], g["in_prompt"][s:e],
            g["gx"][s:e], g["gy"][s:e], g["gz"][s:e],
            sigma_t_ns, rng, vois, outline, args.matching)
        if out is None:
            continue
        bs, bt, br = out
        for i in acc_s:
            acc_s[i].extend(bs[i]); acc_t[i].extend(bt[i]); acc_r[i].extend(br[i])
        if k % 32 == 0:
            print(f"  shard {k}/{len(bounds)}", flush=True)
    f.close()

    rows = []
    for i in [-1] + list(range(len(vois))):
        ts = np.concatenate(acc_s[i]) if acc_s[i] else np.array([], np.float32)
        tt = np.concatenate(acc_t[i]) if acc_t[i] else np.array([], bool)
        tr = np.concatenate(acc_r[i]) if acc_r[i] else np.array([], np.float32)
        n_total, n_sub, n_true_sig, n_true_acc = counts_in_fit_window(ts, tt, tr)
        over = n_sub - n_true_acc
        phi = over / n_true_sig if n_true_sig > 0 else float("nan")
        # counting error on the difference of two independent Poisson counts
        d_over = float(np.sqrt(max(n_sub, 0) + max(n_true_acc, 0)))
        name = "bg" if i == -1 else vois[i]["name"]
        rows.append(dict(voi=name, tau_cfg_ns=(None if i == -1 else vois[i].get("tau_cfg_ns")),
                         n_total=n_total, n_sub=n_sub, n_true_acc=n_true_acc,
                         n_true_sig=n_true_sig, over_subtraction=over,
                         over_sigma=d_over, phi=phi,
                         over_over_sigma=(over / d_over if d_over > 0 else float("nan")),
                         true_acc_frac=(n_true_acc / n_total if n_total else float("nan"))))
        print(f"  {name:>20}  n_tot={n_total:>8,} n_sub={n_sub:>8,} n_true_acc={n_true_acc:>8,} "
              f"n_true_sig={n_true_sig:>8,}  over={over:>+8,} ({over/d_over if d_over>0 else 0:+.1f} sd)"
              f"  phi={phi:+.4f}", flush=True)

    json.dump(dict(label=args.label, matching=args.matching, phantom=args.phantom,
                   singles=str(args.singles), tau_min_ns=EX.TAU_MIN_NS,
                   fit_hi_ns=EX.W_PROMPT_NS, delay_ns=EX.DELAY_NS, rows=rows),
              open(args.out, "w"), indent=1)
    print(f"[{args.label}/{args.matching}] wrote {args.out}")


if __name__ == "__main__":
    main()
