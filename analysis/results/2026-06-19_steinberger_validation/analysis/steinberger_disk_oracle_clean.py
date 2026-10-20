"""Steinberger disk oracle: in-disk-only re-extraction.

The default oracle (steinberger_extract.py) uses the TRUTH SOURCE position
(sx, sy, sz) to truth-correct τ. For thin disk phantoms in air, ⁴⁴Sc β⁺
with σ_β+ ≈ 2 mm in dense glass routinely escapes the 5 mm-thick disk
through the edge or face and annihilates in air or LSO. The source
position is still inside the disk, but the actual annihilation is far
away, and the truth-source-corrected formula

    τ_oracle = delay + (d_annih_to_p1 − d_src_to_p1) / c

picks up a geometric mishandling of order (d_annih − d_src)/c, which
biases the oracle cumulant +0.3 ns.

This script removes those escape events by cutting on the **truth
annihilation position**. In the oracle (no σ_T smear) the TOF-
reconstructed P_ann_truth = the true annihilation position exactly
(by the LOR-midpoint + Δt/c closure). So requiring P_ann_truth ∈ disk
volume IS a truth material cut, same physical meaning as a
`/PLI_LIFETIME: <material>` tag, computed from the existing Singles
HDF5 with no re-simulation.

Output (per disk):
  - N_pairs_total, N_in_disk_pairs, escaped_fraction
  - oracle τ̂ on surviving (in-disk) events
  - bias vs tau_cfg and vs Steinberger reference
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import h5py
import numpy as np

C_MM_NS = 299.792458
TAU_MIN_NS = 1.0

DISKS = [
    {"name": "quartz",        "cx": -50.0, "cy": 0.0, "cz": 0.0,
     "R_mm": 12.7, "halfz_mm": 2.5, "tau_cfg": 1.56,
     "tau_steinberger_label": "1.56 (Takyu) / 1.60 (clinical-isotope cluster)"},
    {"name": "polycarbonate", "cx": +50.0, "cy": 0.0, "cz": 0.0,
     "R_mm":  5.0, "halfz_mm": 2.5, "tau_cfg": 2.07,
     "tau_steinberger_label": "2.07 (Takyu) / 1.927 (²²Na, laminate-caveated)"},
]


def find_shard_boundaries(shard_id):
    if len(shard_id) == 0: return []
    b = np.flatnonzero(np.diff(shard_id) != 0) + 1
    s = np.concatenate(([0], b))
    e = np.concatenate((b, [len(shard_id)]))
    return [(int(shard_id[i]), int(i), int(j)) for i, j in zip(s, e)]


def tof_p_ann(p1, p2, t1, t2):
    mp = 0.5 * (p1 + p2)
    lor = p2 - p1
    norm = np.linalg.norm(lor, axis=-1, keepdims=True)
    direction = lor / np.maximum(norm, 1e-9)
    offset = 0.5 * C_MM_NS * (t1 - t2)[:, None]
    return mp + offset * direction


def per_event_tau_p(p1, p2, pp, t1, t2, tp, p_ann):
    d1 = np.linalg.norm(p1 - p_ann, axis=-1)
    dp = np.linalg.norm(pp - p_ann, axis=-1)
    return (t1 - d1 / C_MM_NS) - (tp - dp / C_MM_NS)


def in_disk(pos, disk):
    dx = pos[:, 0] - disk["cx"]
    dy = pos[:, 1] - disk["cy"]
    dz = pos[:, 2] - disk["cz"]
    return (dx*dx + dy*dy <= disk["R_mm"]**2) & (np.abs(dz) <= disk["halfz_mm"])


def oracle_per_shard(eid, t, gx, gy, gz, sx, sy, sz, in_511, in_prompt, disks):
    """Per shard: oracle-pair by eventID, compute P_ann_truth, per-disk
    return (τ_kept, τ_all) plus counts so we can report escaped fraction."""
    n = len(eid)
    if n == 0: return None
    order = np.argsort(eid, kind="stable")
    eid_s = eid[order]; t_s = t[order]
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

    p1 = np.stack([gx_s[a1], gy_s[a1], gz_s[a1]], axis=-1)
    p2 = np.stack([gx_s[a2], gy_s[a2], gz_s[a2]], axis=-1)
    pp = np.stack([gx_s[p],  gy_s[p],  gz_s[p]],  axis=-1)
    src = np.stack([sx_s[a1], sy_s[a1], sz_s[a1]], axis=-1)
    t1 = t_s[a1]; t2 = t_s[a2]; tp = t_s[p]

    # TRUTH annihilation position: TOF closure with TRUE times = exact annih
    p_ann_truth = tof_p_ann(p1, p2, t1, t2)

    # ALL-events oracle τ (truth-source-corrected, biased by escape geometry)
    tau_all_src = per_event_tau_p(p1, p2, pp, t1, t2, tp, src)
    # In-disk oracle τ (still truth-source-corrected: in-disk means
    # src ≈ annih within positron range so this gives τ ≈ delay)
    tau_all_pann = per_event_tau_p(p1, p2, pp, t1, t2, tp, p_ann_truth)

    # Per disk: route pairs by SOURCE position (parent ion's home disk)
    out = {d["name"]: {"tau_src_all": [], "tau_src_indisk": [],
                       "tau_pann_indisk": [],
                       "n_all": 0, "n_indisk": 0} for d in disks}
    for d in disks:
        sel_src_in = in_disk(src, d)
        sel_annih_in = in_disk(p_ann_truth, d)
        # Pairs where the parent ion was in THIS disk
        n_all_d = int(sel_src_in.sum())
        # Of those, pairs where the annihilation happened IN this disk
        sel_in = sel_src_in & sel_annih_in
        n_in_d = int(sel_in.sum())
        out[d["name"]]["n_all"] += n_all_d
        out[d["name"]]["n_indisk"] += n_in_d
        if sel_src_in.any():
            out[d["name"]]["tau_src_all"].append(tau_all_src[sel_src_in].astype(np.float32))
        if sel_in.any():
            out[d["name"]]["tau_src_indisk"].append(tau_all_src[sel_in].astype(np.float32))
            out[d["name"]]["tau_pann_indisk"].append(tau_all_pann[sel_in].astype(np.float32))
    return out


def cumulant(arr):
    sel = (arr > TAU_MIN_NS) & np.isfinite(arr)
    n = int(sel.sum())
    if n == 0:
        return float("nan"), float("nan"), 0
    tau_hat = float((arr[sel] - TAU_MIN_NS).mean())
    sigma = tau_hat / np.sqrt(max(n, 1))
    return tau_hat, sigma, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--singles", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--max-shards", type=int, default=None)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"[disk-oracle] opening {args.singles}", flush=True)
    f = h5py.File(args.singles, "r"); g = f["singles"]
    sh = g["shard_id"][:]
    bounds = find_shard_boundaries(sh)
    if args.max_shards: bounds = bounds[:args.max_shards]
    print(f"[disk-oracle] {len(bounds)} shards", flush=True)

    acc = {d["name"]: {"tau_src_all": [], "tau_src_indisk": [],
                       "tau_pann_indisk": [],
                       "n_all": 0, "n_indisk": 0} for d in DISKS}
    for sh_id, s, e in bounds:
        eid = g["eventID"][s:e]
        t = g["time_ns"][s:e]
        gx = g["gx"][s:e]; gy = g["gy"][s:e]; gz = g["gz"][s:e]
        sx = g["sx"][s:e]; sy = g["sy"][s:e]; sz = g["sz"][s:e]
        in_5 = g["in_511"][s:e]; in_p = g["in_prompt"][s:e]
        out = oracle_per_shard(eid, t, gx, gy, gz, sx, sy, sz, in_5, in_p, DISKS)
        del eid, t, gx, gy, gz, sx, sy, sz, in_5, in_p
        if out is None: continue
        for k, v in out.items():
            acc[k]["tau_src_all"].extend(v["tau_src_all"])
            acc[k]["tau_src_indisk"].extend(v["tau_src_indisk"])
            acc[k]["tau_pann_indisk"].extend(v["tau_pann_indisk"])
            acc[k]["n_all"]    += v["n_all"]
            acc[k]["n_indisk"] += v["n_indisk"]
    f.close()

    rows = []
    print(f"\n[disk-oracle] per-disk in-disk re-extraction:")
    print(f"{'disk':>16} {'N_pairs':>10} {'N_in_disk':>10} {'esc_frac':>10}  "
          f"{'τ̂(all,src)':>12} {'τ̂(in,src)':>12} {'τ̂(in,P_ann)':>14}  "
          f"{'cfg':>6} {'bias_vs_cfg_ps':>15}")
    for d in DISKS:
        a = acc[d["name"]]
        arr_all = np.concatenate(a["tau_src_all"]) if a["tau_src_all"] else np.array([], dtype=np.float32)
        arr_in_src = np.concatenate(a["tau_src_indisk"]) if a["tau_src_indisk"] else np.array([], dtype=np.float32)
        arr_in_p = np.concatenate(a["tau_pann_indisk"]) if a["tau_pann_indisk"] else np.array([], dtype=np.float32)
        tau_all_src, sig_all_src, n_all_tail   = cumulant(arr_all)
        tau_in_src,  sig_in_src,  n_in_src_tail = cumulant(arr_in_src)
        tau_in_p,    sig_in_p,    n_in_p_tail  = cumulant(arr_in_p)
        escaped = 1.0 - a["n_indisk"] / max(a["n_all"], 1)
        bias_in_p_ps = (tau_in_p - d["tau_cfg"]) * 1000 if np.isfinite(tau_in_p) else float("nan")
        row = dict(
            disk=d["name"], tau_cfg_ns=d["tau_cfg"],
            tau_steinberger_label=d["tau_steinberger_label"],
            N_pairs=a["n_all"], N_in_disk=a["n_indisk"],
            escaped_fraction=escaped,
            oracle_all_src=tau_all_src, oracle_all_src_sigma=sig_all_src,
            oracle_in_disk_src=tau_in_src, oracle_in_disk_src_sigma=sig_in_src,
            oracle_in_disk_pann=tau_in_p,  oracle_in_disk_pann_sigma=sig_in_p,
            N_tail_all=n_all_tail, N_tail_in=n_in_p_tail,
            bias_in_pann_vs_cfg_ps=bias_in_p_ps,
        )
        rows.append(row)
        print(f"{d['name']:>16} {a['n_all']:>10,} {a['n_indisk']:>10,} "
              f"{escaped*100:>9.1f}% "
              f"{tau_all_src:>12.3f} {tau_in_src:>12.3f} {tau_in_p:>14.3f}  "
              f"{d['tau_cfg']:>6.2f} {bias_in_p_ps:>+14.1f}")

    # Save
    import csv
    with open(args.out / "disk_oracle_in_disk.csv", "w") as fp:
        w = csv.DictWriter(fp, fieldnames=rows[0].keys())
        w.writeheader()
        for r in rows: w.writerow(r)
    (args.out / "disk_oracle_in_disk.json").write_text(json.dumps(rows, indent=2,
        default=lambda o: float(o) if isinstance(o, np.floating) else o))
    print(f"\nwrote {args.out}/disk_oracle_in_disk.{{csv,json}}")


if __name__ == "__main__":
    sys.exit(main())
