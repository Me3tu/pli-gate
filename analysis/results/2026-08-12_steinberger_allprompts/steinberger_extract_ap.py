"""Steinberger validation: dual-pipeline τ̂ extractor.

Runs BOTH pipelines on the same compact Singles HDF5:
  (a) ORACLE  : per-shard eventID grouping, truth-position-corrected τ,
                cumulant MLE on τ > τ_min. Geometric VOI assignment by
                source-truth position (sx, sy, sz).
  (b) CLINICAL: M3_tubes pattern. Sliding ±2 ns pair sorter + ±5 ns
                prompt window + +50 ns delayed-window randoms + TOF-
                corrected τ + truncation-aware MLE. Geometric VOI
                assignment by reconstructed annihilation point (p_ann).

VOIs come from a JSON config with shape ∈ {"sphere", "cylinder"}.

Per-VOI output:
  oracle:   {tau_hat_ns, sigma_ns, N_signal}
  clinical: {tau_hat_naive_ns, tau_hat_trunc_ns, sigma_ns, sigma_trunc_ns,
             N_total, N_random, N_signal, detection_sigma_vs_bg}

Compare against tau_cfg_ns (configured input) AND tau_steinberger_ns
(from the paper).

USAGE
-----
  source ../../../scripts/env_setup.sh   # for the venv
  python3 steinberger_extract.py \
      --singles  <converted Singles HDF5> \
      --config   ../analysis/voi_configs.json \
      --phantom  disk_phantom | water_cylinder \
      --out      <output dir> \
      --sigma-t-ps 90.9
"""
from __future__ import annotations
import argparse, csv, json, sys, time
from pathlib import Path
import h5py
import numpy as np

C_MM_NS = 299.792458
TAU_MIN_NS = 1.0

W_PAIR_NS    = 2.0   # clinical
W_PROMPT_NS  = 5.0   # clinical
DELAY_NS     = 50.0  # clinical
TAU_BINS     = np.linspace(-5.0, 8.0, 261)


# =======================================================================
# Shared helpers
# =======================================================================

def find_shard_boundaries(shard_id):
    if len(shard_id) == 0: return []
    b = np.flatnonzero(np.diff(shard_id) != 0) + 1
    s = np.concatenate(([0], b))
    e = np.concatenate((b, [len(shard_id)]))
    return [(int(shard_id[i]), int(i), int(j)) for i, j in zip(s, e)]


def in_phantom_outline(x, y, z, outline):
    if outline["type"] == "cylinder":
        R = outline["R_mm"] + outline.get("buffer_mm", 0)
        H = outline["halfz_mm"] + outline.get("buffer_mm", 0)
        return ((x*x + y*y) <= R*R) & (np.abs(z) <= H)
    elif outline["type"] == "world_box":
        hx = outline["halfx_mm"] + outline.get("buffer_mm", 0)
        hy = outline["halfy_mm"] + outline.get("buffer_mm", 0)
        hz = outline["halfz_mm"] + outline.get("buffer_mm", 0)
        return (np.abs(x) <= hx) & (np.abs(y) <= hy) & (np.abs(z) <= hz)
    else:
        raise ValueError(f"unknown phantom outline type: {outline['type']}")


def voi_masks(pos, vois):
    """Return a list of N boolean masks, one per VOI. Each mask is True for
    triples geometrically inside that VOI. VOIs may overlap (nested-cylinder
    Steinberger pattern). Each triple is independently counted in EACH VOI it
    satisfies. Bg mask = NOT in any VOI.
    """
    masks = []
    for v in vois:
        dx = pos[:, 0] - v["cx_mm"]
        dy = pos[:, 1] - v["cy_mm"]
        dz = pos[:, 2] - v["cz_mm"]
        if v["shape"] == "sphere":
            d2 = dx*dx + dy*dy + dz*dz
            sel = d2 <= v["r_mm"] * v["r_mm"]
        elif v["shape"] == "cylinder":
            d_xy2 = dx*dx + dy*dy
            sel = (d_xy2 <= v["r_mm"]*v["r_mm"]) & (np.abs(dz) <= v["half_z_mm"])
        else:
            raise ValueError(f"unknown VOI shape: {v['shape']}")
        masks.append(sel)
    return masks


def voi_bg_mask(pos, vois):
    """bg = NOT in any VOI."""
    in_any = np.zeros(len(pos), dtype=bool)
    for m in voi_masks(pos, vois):
        in_any |= m
    return ~in_any


def lor_perp(p_ann, p1, p2):
    lor = p2 - p1
    norm = np.linalg.norm(lor, axis=-1, keepdims=True)
    direction = lor / np.maximum(norm, 1e-9)
    v = p_ann - p1
    proj = (v * direction).sum(axis=-1, keepdims=True) * direction
    return np.linalg.norm(v - proj, axis=-1)


def tof_p_ann(p1, p2, t1, t2):
    mp = 0.5 * (p1 + p2)
    lor = p2 - p1
    norm = np.linalg.norm(lor, axis=-1, keepdims=True)
    direction = lor / np.maximum(norm, 1e-9)
    offset = 0.5 * C_MM_NS * (t1 - t2)[:, None]
    return mp + offset * direction


def per_event_tau(p1, p2, pp, t1, t2, tp, p_ann):
    d1 = np.linalg.norm(p1 - p_ann, axis=-1)
    dp = np.linalg.norm(pp - p_ann, axis=-1)
    return (t1 - d1 / C_MM_NS) - (tp - dp / C_MM_NS)


def closest_prompt_indices(t_centers, t_prompt_sorted, half_ns):
    lo = np.searchsorted(t_prompt_sorted, t_centers - half_ns, side="left")
    hi = np.searchsorted(t_prompt_sorted, t_centers + half_ns, side="right")
    has = hi > lo
    cand_lo = lo[has]; cand_hi = np.maximum(hi[has] - 1, cand_lo)
    d_lo = np.abs(t_prompt_sorted[cand_lo] - t_centers[has])
    d_hi = np.abs(t_prompt_sorted[cand_hi] - t_centers[has])
    pick = np.where(d_lo <= d_hi, cand_lo, cand_hi)
    return has, pick


def all_prompt_indices(t_centers, t_prompt_sorted, half_ns):
    """S1 port from Phase 5 extract_sim1_scatter_ewin.py: every prompt in
    the window forms a triple (accidentals exactly additive, so the
    delayed-window randoms estimate is unbiased)."""
    lo = np.searchsorted(t_prompt_sorted, t_centers - half_ns, side="left")
    hi = np.searchsorted(t_prompt_sorted, t_centers + half_ns, side="right")
    counts = hi - lo
    pair_idx = np.repeat(np.arange(len(t_centers)), counts)
    pick = (np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)])
            if counts.sum() else np.array([], dtype=np.int64))
    return pair_idx, pick


def truncation_correct(tau_obs, W):
    from scipy.optimize import brentq
    if not np.isfinite(tau_obs) or tau_obs <= 0 or W <= 0: return float("nan")
    def f(lam):
        u = W/lam
        if u > 50: return lam - tau_obs
        return lam - W*np.exp(-u)/(1-np.exp(-u)) - tau_obs
    try:    return float(brentq(f, tau_obs*1.001, 50.0, xtol=1e-6))
    except: return float("nan")


# =======================================================================
# Oracle pipeline (eventID grouping, truth-position τ)
# =======================================================================

def oracle_per_shard(eid, time_ns, gx, gy, gz, sx, sy, sz, in_511, in_prompt,
                    vois):
    """Returns by_voi dict: {voi_idx: [tau_array]}."""
    n = len(eid)
    if n == 0: return None
    order = np.argsort(eid, kind="stable")
    eid_s = eid[order]; t_s = time_ns[order]
    gx_s, gy_s, gz_s = gx[order], gy[order], gz[order]
    sx_s, sy_s, sz_s = sx[order], sy[order], sz[order]
    in_511_s = in_511[order].astype(bool)
    in_prompt_s = in_prompt[order].astype(bool)
    is_new = np.empty(n, dtype=bool); is_new[0] = True
    is_new[1:] = eid_s[1:] != eid_s[:-1]
    group_idx = np.cumsum(is_new) - 1
    n_groups = int(group_idx[-1]) + 1
    n_511_g = np.bincount(group_idx[in_511_s], minlength=n_groups)
    n_prompt_g = np.bincount(group_idx[in_prompt_s], minlength=n_groups)
    valid_g = (n_511_g == 2) & (n_prompt_g >= 1)
    n_valid = int(valid_g.sum())
    if n_valid == 0: return None
    valid_511 = valid_g[group_idx] & in_511_s
    valid_prompt = valid_g[group_idx] & in_prompt_s
    idx_511_v = np.flatnonzero(valid_511)
    if len(idx_511_v) != 2 * n_valid: return None
    pair_order = np.argsort(group_idx[idx_511_v], kind="stable")
    paired = idx_511_v[pair_order].reshape(-1, 2)
    a1, a2 = paired[:, 0], paired[:, 1]
    t_mean = 0.5 * (t_s[a1] + t_s[a2])
    valid_row = np.full(n_groups, -1, dtype=np.int64)
    valid_row[np.flatnonzero(valid_g)] = np.arange(n_valid)
    idx_p_v = np.flatnonzero(valid_prompt)
    row_p = valid_row[group_idx[idx_p_v]]
    dt_p = np.abs(t_s[idx_p_v] - t_mean[row_p])
    so = np.lexsort((dt_p, row_p)); rs = row_p[so]
    first = np.empty(len(rs), dtype=bool); first[0] = True
    first[1:] = rs[1:] != rs[:-1]
    pick = so[first]
    p = np.empty(n_valid, dtype=np.int64); p[rs[first]] = idx_p_v[pick]

    # Truth-position-corrected τ (uses source position, not TOF)
    src = np.stack([sx_s[a1], sy_s[a1], sz_s[a1]], axis=-1)
    p1 = np.stack([gx_s[a1], gy_s[a1], gz_s[a1]], axis=-1)
    p2 = np.stack([gx_s[a2], gy_s[a2], gz_s[a2]], axis=-1)
    pp = np.stack([gx_s[p],  gy_s[p],  gz_s[p]],  axis=-1)
    tau = per_event_tau(p1, p2, pp, t_s[a1], t_s[a2], t_s[p], src)
    # VOI assignment by source truth position, INDEPENDENT per VOI
    # (each triple may be counted in multiple nested VOIs)
    masks = voi_masks(src, vois)
    by_voi = {i: [] for i in range(-1, len(vois))}
    for i, m in enumerate(masks):
        if m.any():
            by_voi[i].append(tau[m].astype(np.float32))
    bg_mask = voi_bg_mask(src, vois)
    if bg_mask.any():
        by_voi[-1].append(tau[bg_mask].astype(np.float32))
    return by_voi


def oracle_run(singles_h5, vois, sigma_t_ps=0.0, max_shards=None):
    f = h5py.File(singles_h5, "r"); g = f["singles"]
    sh = g["shard_id"][:]
    bounds = find_shard_boundaries(sh)
    if max_shards: bounds = bounds[:max_shards]
    acc = {i: [] for i in range(-1, len(vois))}
    # For the ORACLE, sigma_t_ps is OPTIONAL. If 0, no smear (matches
    # M4 truth-time recovery). If >0, applies σ_T smear to study CRT.
    rng = np.random.default_rng(int(sigma_t_ps * 1000) or 1)
    sigma_t_ns = sigma_t_ps / 1000.0
    for sh_id, s, e in bounds:
        eid = g["eventID"][s:e]
        t_true = g["time_ns"][s:e]
        if sigma_t_ps > 0:
            t_true = t_true + rng.normal(0.0, sigma_t_ns, size=len(t_true))
        gx = g["gx"][s:e]; gy = g["gy"][s:e]; gz = g["gz"][s:e]
        sx = g["sx"][s:e]; sy = g["sy"][s:e]; sz = g["sz"][s:e]
        in_5 = g["in_511"][s:e]; in_p = g["in_prompt"][s:e]
        out = oracle_per_shard(eid, t_true, gx, gy, gz, sx, sy, sz, in_5, in_p, vois)
        del eid, t_true, gx, gy, gz, sx, sy, sz, in_5, in_p
        if out:
            for k, lst in out.items(): acc[k].extend(lst)
    f.close()
    # Cumulant MLE per VOI
    rows = []
    for k in [-1] + list(range(len(vois))):
        arr = np.concatenate(acc[k]) if acc[k] else np.array([], dtype=np.float32)
        sel = (arr > TAU_MIN_NS) & np.isfinite(arr)
        if not sel.any():
            rows.append(dict(voi=("bg" if k==-1 else vois[k]["name"]),
                              tau_hat_ns=float("nan"), sigma_ns=float("nan"),
                              n_signal=0))
            continue
        tau_hat = float((arr[sel] - TAU_MIN_NS).mean())
        n = int(sel.sum())
        sigma = tau_hat / np.sqrt(max(n, 1))
        rows.append(dict(voi=("bg" if k==-1 else vois[k]["name"]),
                          tau_hat_ns=tau_hat, sigma_ns=sigma, n_signal=n))
    return rows


# =======================================================================
# Clinical pipeline (sliding pair + delayed-window randoms + trunc MLE)
# =======================================================================

def clinical_per_shard(t_true, in_511, in_prompt, gx, gy, gz, sigma_t_ns, rng,
                       vois, phantom_outline, matching="closest"):
    n = len(t_true)
    if n == 0: return None
    t_smear = t_true + rng.normal(0.0, sigma_t_ns, size=n)
    in_511_b = in_511.astype(bool); in_prompt_b = in_prompt.astype(bool)
    idx_511 = np.flatnonzero(in_511_b); idx_prompt = np.flatnonzero(in_prompt_b)
    if len(idx_511) < 2 or len(idx_prompt) < 1: return None
    t_511 = t_smear[idx_511]; t_prompt = t_smear[idx_prompt]
    o5 = np.argsort(t_511); op = np.argsort(t_prompt)
    idx_511 = idx_511[o5]; t_511 = t_511[o5]
    idx_prompt = idx_prompt[op]; t_prompt = t_prompt[op]
    dt = np.diff(t_511); pair_ok = dt <= W_PAIR_NS
    if not pair_ok.any(): return None
    pa = idx_511[:-1][pair_ok]; pb = idx_511[1:][pair_ok]
    p1 = np.stack([gx[pa], gy[pa], gz[pa]], -1)
    p2 = np.stack([gx[pb], gy[pb], gz[pb]], -1)
    p_ann = tof_p_ann(p1, p2, t_smear[pa], t_smear[pb])
    perp = lor_perp(p_ann, p1, p2)
    geom = (perp <= phantom_outline.get("lor_tol_mm", 10.0)) & in_phantom_outline(
        p_ann[:,0], p_ann[:,1], p_ann[:,2], phantom_outline)
    if not geom.any(): return None
    pa = pa[geom]; pb = pb[geom]; p_ann = p_ann[geom]
    t_pair_mean = 0.5*(t_smear[pa] + t_smear[pb])

    if matching == "all":
        sel_s, pick_s = all_prompt_indices(t_pair_mean, t_prompt, W_PROMPT_NS)
        sel_r, pick_r = all_prompt_indices(t_pair_mean + DELAY_NS, t_prompt, W_PROMPT_NS)
    else:
        has_s, pick_s = closest_prompt_indices(t_pair_mean, t_prompt, W_PROMPT_NS)
        has_r, pick_r = closest_prompt_indices(t_pair_mean + DELAY_NS, t_prompt, W_PROMPT_NS)
        sel_s = np.flatnonzero(has_s); sel_r = np.flatnonzero(has_r)

    def build_tau(sel, pick):
        if len(sel) == 0: return None, None
        pa_h = pa[sel]; pb_h = pb[sel]; p_ann_h = p_ann[sel]
        ip = idx_prompt[pick]
        p1_h = np.stack([gx[pa_h], gy[pa_h], gz[pa_h]], -1)
        p2_h = np.stack([gx[pb_h], gy[pb_h], gz[pb_h]], -1)
        pp_h = np.stack([gx[ip], gy[ip], gz[ip]], -1)
        tau = per_event_tau(p1_h, p2_h, pp_h,
                            t_smear[pa_h], t_smear[pb_h], t_smear[ip], p_ann_h)
        return tau, p_ann_h
    tau_s, p_ann_s = build_tau(sel_s, pick_s)
    tau_r, p_ann_r = build_tau(sel_r, pick_r)
    if tau_r is not None: tau_r = tau_r + DELAY_NS

    by_sig = {i: [] for i in range(-1, len(vois))}
    by_rand = {i: [] for i in range(-1, len(vois))}
    if tau_s is not None:
        masks_s = voi_masks(p_ann_s, vois)
        for i, m in enumerate(masks_s):
            if m.any(): by_sig[i].append(tau_s[m].astype(np.float32))
        bg_s = voi_bg_mask(p_ann_s, vois)
        if bg_s.any(): by_sig[-1].append(tau_s[bg_s].astype(np.float32))
    if tau_r is not None:
        masks_r = voi_masks(p_ann_r, vois)
        for i, m in enumerate(masks_r):
            if m.any(): by_rand[i].append(tau_r[m].astype(np.float32))
        bg_r = voi_bg_mask(p_ann_r, vois)
        if bg_r.any(): by_rand[-1].append(tau_r[bg_r].astype(np.float32))
    return by_sig, by_rand


def subtract_and_mle(tau_total, tau_random):
    h_t, _ = np.histogram(tau_total, bins=TAU_BINS)
    h_r, _ = np.histogram(tau_random, bins=TAU_BINS)
    bc = 0.5*(TAU_BINS[:-1] + TAU_BINS[1:])
    sel = (bc > TAU_MIN_NS) & (bc < W_PROMPT_NS)
    h_t_t = h_t[sel].astype(np.float64); h_r_t = h_r[sel].astype(np.float64)
    bc_t = bc[sel]
    h_s = h_t_t - h_r_t
    n_s = float(h_s.sum())
    if n_s <= 0:
        return dict(tau_naive=float("nan"), tau_trunc=float("nan"),
                    sigma_naive=float("nan"), sigma_trunc=float("nan"),
                    n_total=int(h_t_t.sum()), n_random=int(h_r_t.sum()),
                    n_signal=n_s)
    tau_obs = float((bc_t * h_s).sum() / n_s) - TAU_MIN_NS
    W = W_PROMPT_NS - TAU_MIN_NS
    tau_trunc = truncation_correct(tau_obs, W)
    n_var = h_t_t.sum() + h_r_t.sum()
    sig_naive = float(abs(tau_obs)*np.sqrt(n_var)/max(n_s, 1))
    if np.isfinite(tau_trunc):
        eps = 1e-3
        lam_p = truncation_correct(tau_obs+eps, W)
        sig_trunc = float(abs((lam_p - tau_trunc)/eps)*sig_naive) if np.isfinite(lam_p) else float("nan")
    else:
        sig_trunc = float("nan")
    return dict(tau_naive=tau_obs, tau_trunc=tau_trunc,
                sigma_naive=sig_naive, sigma_trunc=sig_trunc,
                n_total=int(h_t_t.sum()), n_random=int(h_r_t.sum()),
                n_signal=n_s)


def clinical_run(singles_h5, vois, phantom_outline, sigma_t_ps, matching="closest",
                 max_shards=None):
    f = h5py.File(singles_h5, "r"); g = f["singles"]
    sh = g["shard_id"][:]; bounds = find_shard_boundaries(sh)
    if max_shards: bounds = bounds[:max_shards]
    rng = np.random.default_rng(int(sigma_t_ps * 1000) or 1)
    sigma_t_ns = sigma_t_ps / 1000.0
    acc_sig  = {i: [] for i in range(-1, len(vois))}
    acc_rand = {i: [] for i in range(-1, len(vois))}
    for sh_id, s, e in bounds:
        t_true = g["time_ns"][s:e]
        in_5 = g["in_511"][s:e]; in_p = g["in_prompt"][s:e]
        gx = g["gx"][s:e]; gy = g["gy"][s:e]; gz = g["gz"][s:e]
        out = clinical_per_shard(t_true, in_5, in_p, gx, gy, gz,
                                  sigma_t_ns, rng, vois, phantom_outline,
                                  matching=matching)
        del t_true, in_5, in_p, gx, gy, gz
        if out is None: continue
        by_sig, by_rand = out
        for k, lst in by_sig.items():  acc_sig[k].extend(lst)
        for k, lst in by_rand.items(): acc_rand[k].extend(lst)
    f.close()
    rows = []
    for k in [-1] + list(range(len(vois))):
        sig = np.concatenate(acc_sig[k]) if acc_sig[k] else np.array([], dtype=np.float32)
        rnd = np.concatenate(acc_rand[k]) if acc_rand[k] else np.array([], dtype=np.float32)
        res = subtract_and_mle(sig, rnd)
        rows.append(dict(voi=("bg" if k==-1 else vois[k]["name"]),
                          **res))
    return rows


# =======================================================================
# Main
# =======================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--singles", required=True, type=Path)
    ap.add_argument("--config",  required=True, type=Path)
    # The choice list used to be hard-coded to ["disk_phantom",
    # "water_cylinder"], which meant a new phantom could not be extracted
    # without editing this file. It is now validated against the config that
    # was actually passed, which is both more permissive and stricter: a typo
    # is still rejected, and the error names the keys the file really has.
    ap.add_argument("--phantom", required=True)
    ap.add_argument("--out",     required=True, type=Path)
    ap.add_argument("--sigma-t-ps", type=float, default=90.9)
    ap.add_argument("--matching", choices=["closest", "all"], default="closest")
    ap.add_argument("--max-shards", type=int, default=None)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    all_cfg = json.loads(args.config.read_text())
    if args.phantom not in all_cfg:
        raise SystemExit(
            f"--phantom {args.phantom!r} is not a key of {args.config}. "
            f"Available: {sorted(all_cfg)}")
    cfg = all_cfg[args.phantom]
    vois = cfg["vois"]; outline = cfg["phantom_outline"]

    t0 = time.time()
    print(f"[steinberger] {args.phantom}  σ_T={args.sigma_t_ps} ps", flush=True)
    print(f"[steinberger] VOIs ({len(vois)}):")
    for v in vois:
        print(f"   {v['name']:>16}  shape={v['shape']:<8} "
              f"@ ({v['cx_mm']:+.1f}, {v['cy_mm']:+.1f}, {v['cz_mm']:+.1f}) mm  "
              f"tau_cfg={v.get('tau_cfg_ns')}  tau_steinberger={v.get('tau_steinberger_ns', v.get('tau_steinberger_ns_range'))}")

    print(f"\n[steinberger] running ORACLE pipeline (truth-position-corrected τ, no σ_T smear)")
    oracle_rows = oracle_run(args.singles, vois, sigma_t_ps=0.0,
                              max_shards=args.max_shards)
    print(f"[steinberger] running CLINICAL pipeline (σ_T={args.sigma_t_ps} ps re-smear)")
    clinical_rows = clinical_run(args.singles, vois, outline,
                                  args.sigma_t_ps, matching=args.matching,
                                  max_shards=args.max_shards)

    # Merge into a single table
    oracle_by_name = {r["voi"]: r for r in oracle_rows}
    clinical_by_name = {r["voi"]: r for r in clinical_rows}
    table = []
    for v in vois + [{"name":"bg"}]:
        name = v["name"]
        o = oracle_by_name.get(name, {})
        c = clinical_by_name.get(name, {})
        st = v.get("tau_steinberger_ns") or v.get("tau_steinberger_ns_range")
        row = dict(
            voi=name,
            tau_cfg_ns=v.get("tau_cfg_ns"),
            tau_steinberger=st,
            oracle_tau_hat_ns=o.get("tau_hat_ns"),
            oracle_sigma_ns=o.get("sigma_ns"),
            oracle_n_signal=o.get("n_signal"),
            clinical_tau_naive_ns=c.get("tau_naive"),
            clinical_tau_trunc_ns=c.get("tau_trunc"),
            clinical_sigma_naive_ns=c.get("sigma_naive"),
            clinical_sigma_trunc_ns=c.get("sigma_trunc"),
            clinical_n_total=c.get("n_total"),
            clinical_n_random=c.get("n_random"),
            clinical_n_signal=c.get("n_signal"),
        )
        table.append(row)

    # Print
    print(f"\n[steinberger] τ̂ comparison table:")
    print(f"{'VOI':>16} {'cfg':>5} {'Stein.':>10} {'oracle':>8} {'σ_o':>7} "
          f"{'clin_trunc':>11} {'σ_c':>7} {'N_clin_sig':>11}")
    for r in table:
        cfg_s = f"{r['tau_cfg_ns']:.2f}" if r['tau_cfg_ns'] is not None else "    -"
        st_s = str(r['tau_steinberger']) if r['tau_steinberger'] else "       -"
        o_s = f"{r['oracle_tau_hat_ns']:.3f}" if r['oracle_tau_hat_ns'] is not None and np.isfinite(r['oracle_tau_hat_ns']) else "      -"
        oσ_s = f"{r['oracle_sigma_ns']:.3f}" if r['oracle_sigma_ns'] is not None and np.isfinite(r['oracle_sigma_ns']) else "    -"
        c_s = f"{r['clinical_tau_trunc_ns']:.3f}" if r['clinical_tau_trunc_ns'] is not None and np.isfinite(r['clinical_tau_trunc_ns']) else "         -"
        cσ_s = f"{r['clinical_sigma_trunc_ns']:.3f}" if r['clinical_sigma_trunc_ns'] is not None and np.isfinite(r['clinical_sigma_trunc_ns']) else "    -"
        n_s = f"{r['clinical_n_signal']:.0f}" if r['clinical_n_signal'] is not None else "       -"
        print(f"{r['voi']:>16} {cfg_s:>5} {st_s:>10} {o_s:>8} {oσ_s:>7} {c_s:>11} {cσ_s:>7} {n_s:>11}")

    # Save
    with open(args.out / "per_voi.csv", "w") as fp:
        if table:
            w = csv.DictWriter(fp, fieldnames=table[0].keys())
            w.writeheader()
            for r in table:
                # JSON-serialize ranges + None
                rr = {k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()}
                w.writerow(rr)
    (args.out / "state.json").write_text(json.dumps(dict(
        phantom=args.phantom, sigma_t_ps=args.sigma_t_ps,
        config=cfg, oracle=oracle_rows, clinical=clinical_rows,
        wall_clock_s=float(time.time() - t0),
    ), default=lambda o: float(o) if isinstance(o, np.floating) else o, indent=2))
    print(f"\n[steinberger] wrote {args.out}/per_voi.csv  state.json")


if __name__ == "__main__":
    sys.exit(main())
