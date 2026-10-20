"""Phase 3s-2 — Hypoxia per-insert MLE tau.

Same MLE pattern as Phase 2b analyze_per_tube_tau.py but with the 5-insert
hypoxia layout (spheres at x=-60,-30,0,+30,+60 mm, all r=15 mm). Each
insert is a different TumorTissue_pO2_N material with a slightly different
tau_TE per the materials JSON.
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
import h5py, numpy as np

REPO = Path(os.environ.get("PLI_REPO") or Path(__file__).resolve().parents[4]).resolve()
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from reconstruction.correct_travel_time import correct as travel_time_correct  # noqa

INSERTS = [
    ("TumorTissue_pO2_2",  -60.0, 0.0, 0.0, 15.0),  # name, cx, cy, cz, r_mm
    ("TumorTissue_pO2_6",  -30.0, 0.0, 0.0, 15.0),
    ("TumorTissue_pO2_10",   0.0, 0.0, 0.0, 15.0),
    ("TumorTissue_pO2_20",  30.0, 0.0, 0.0, 15.0),
    ("TumorTissue_pO2_40",  60.0, 0.0, 0.0, 15.0),
]
ASSIGN_R_MM = 16.0  # 1 mm guard on top of geometric radius


def assign(src_x, src_y, src_z, inserts):
    a = np.full(len(src_x), -1, dtype=np.int8)
    for i, (_, cx, cy, cz, r) in enumerate(inserts):
        d2 = (src_x - cx)**2 + (src_y - cy)**2 + (src_z - cz)**2
        a[d2 <= ASSIGN_R_MM**2] = i
    return a


def mle_tau(tau, tau_min):
    s = tau[(tau > tau_min) & np.isfinite(tau)]
    if s.size < 20:
        return float("nan"), float("nan"), int(s.size)
    t = float(s.mean()) - tau_min
    return t, t / np.sqrt(s.size), int(s.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--detector", required=True, type=Path)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--tau-min", type=float, default=1.0)
    ap.add_argument("--materials-json", type=Path,
                    default=REPO / "materials" / "pli_materials.json")
    ap.add_argument("--out-md", required=True, type=Path)
    ap.add_argument("--out-csv", required=True, type=Path)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    print(f"[hypoxia] travel-time correction")
    tau, _ = travel_time_correct(args.detector)
    print(f"[hypoxia]   N triples={len(tau):,}  median={np.median(tau):.3f} ns  frac>0={(tau>0).mean():.3f}")
    with h5py.File(args.truth, "r") as ft:
        sx = ft["truth"]["source_x"][:]
        sy = ft["truth"]["source_y"][:]
        sz = ft["truth"]["source_z"][:]
    a = assign(sx, sy, sz, INSERTS)
    print(f"[hypoxia] assignment counts (assign_r={ASSIGN_R_MM} mm):")
    for i, (name, cx, cy, cz, r) in enumerate(INSERTS):
        n = int((a == i).sum())
        print(f"     {name:25s} ({cx:+.0f},{cy:+.0f},{cz:+.0f}) r={r}: {n:,}")
    n_un = int((a == -1).sum())
    print(f"     unassigned (water bg): {n_un:,}")

    with open(args.materials_json) as f:
        mats = {m["name"]: m for m in json.load(f)["materials"]}
    rows = []
    print(f"[hypoxia] per-insert MLE (tau_min={args.tau_min}):")
    for i, (name, cx, cy, cz, r) in enumerate(INSERTS):
        m = mats.get(name, {})
        cfg = m.get("tau_TE_ns", float("nan"))
        pub = m.get("tau_pub_ns", float("nan"))
        i3 = m.get("I3", float("nan"))
        idx = (a == i)
        rec, sig, n = mle_tau(tau[idx], args.tau_min)
        resid = (rec - cfg) / cfg * 100 if np.isfinite(rec) and np.isfinite(cfg) and cfg else float("nan")
        n_total = int(idx.sum())
        print(f"     {name:25s}  cfg={cfg:.4f}  rec={rec:.4f}+/-{sig:.4f}  "
              f"resid={resid:+.3f}%  N={n:,}/{n_total:,}")
        rows.append(dict(insert=name, cx=cx, cy=cy, cz=cz, r=r,
                         tau_pub_lit=pub, I3=i3, tau_TE_configured=cfg,
                         tau_recovered=rec, sigma=sig, residual_pct=resid,
                         n_used=n, n_total=n_total))

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    keys = ["insert","cx","cy","cz","r","tau_pub_lit","I3","tau_TE_configured",
            "tau_recovered","sigma","residual_pct","n_used","n_total"]
    with open(args.out_csv, "w") as f:
        f.write(",".join(keys) + "\n")
        for r_ in rows:
            f.write(",".join(str(r_[k]) for k in keys) + "\n")

    # Summary
    valid = [r_ for r_ in rows if np.isfinite(r_["residual_pct"])]
    mean_abs = float(np.mean([abs(r_["residual_pct"]) for r_ in valid])) if valid else float("nan")
    fits = np.array([r_["tau_recovered"] for r_ in rows if np.isfinite(r_["tau_recovered"])])
    spread_ps = float((fits.max() - fits.min()) * 1000) if len(fits) >= 2 else float("nan")

    with open(args.out_md, "w") as f:
        f.write(f"# Hypoxia per-insert τ̂ — {args.label}\n\n")
        f.write(f"Source: `{args.detector.name}`  τ_min={args.tau_min} ns  assign_r={ASSIGN_R_MM} mm\n\n")
        f.write("| insert | tau_pub | I3 | tau_TE_cfg | tau_recovered | sigma | resid % | N |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for r_ in rows:
            f.write(f"| `{r_['insert']}` | {r_['tau_pub_lit']:.4f} | {r_['I3']:.3f} | "
                    f"{r_['tau_TE_configured']:.4f} | {r_['tau_recovered']:.4f} | {r_['sigma']:.4f} | "
                    f"{r_['residual_pct']:+.2f} | {r_['n_used']:,} |\n")
        f.write(f"\n**mean |residual| across {len(valid)} inserts: {mean_abs:.3f} %**\n")
        f.write(f"\n**fit τ̂ spread (max - min): {spread_ps:.2f} ps**\n")
    print(f"[hypoxia] mean |residual| = {mean_abs:.3f} %; spread = {spread_ps:.2f} ps")
    print(f"[hypoxia] wrote {args.out_csv} {args.out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
