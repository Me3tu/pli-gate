#!/usr/bin/env python3
"""MC15 follow-up: is the SpineBone bias present at the level of true annihilations?

Truth oracle (post-correction VIP-Man slice, 10.2 M annihilations: material at the
annihilation site, true delay, true annihilation position). Same estimator as Table 4
(mean(tau-1 | tau>1)), same tissue labelling (modal 4 mm atlas), no detector, no TOF.
  E0  by annihilation MATERIAL                     -> configured lifetime (sanity)
  E1  by atlas label of the 4 mm cell containing the annihilation
  E2  8 mm voxels, modal atlas label, count-weighted (the Table-4 aggregation)
Also: composition of annihilations inside atlas-SpineBone cells, with each material's
share of the tau>1 ns tail.
"""
import gzip, json, os
from collections import Counter
from pathlib import Path
import numpy as np

RES = Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1])
TL = np.load(RES / "simple/run_B_4mm_10iter_geosens_norm/tissue_label_3d.npy", allow_pickle=True)
EXT0 = np.array([-280.0, -140.0, -40.0])
mat, d, pos = [], [], []
if not Path(Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-13_overnight/V2_migration/oracle_records.txt.gz").exists(): raise SystemExit(f"{Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-13_overnight/V2_migration/oracle_records.txt.gz"}: raw data not included in the release, see analysis/README.md")
with gzip.open(Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-08-13_overnight/V2_migration/oracle_records.txt.gz", "rt") as f:
    for line in f:
        p = line.rstrip().split(":")
        if len(p) < 7 or p[0] != "PLI_LIFETIME": continue
        d.append(float(p[1])); mat.append(p[2]); pos.append((float(p[3]), float(p[4]), float(p[5])))
d = np.array(d); mat = np.array(mat); pos = np.array(pos)
print(f"records {len(d):,}")
CFG = {"Adipose": 2.55, "Blood": 1.80, "Body": 2.00, "Cartilage": 2.10, "Heart": 2.00, "Lung": 1.85, "Muscle": 2.03, "SpineBone": 1.60}
est = lambda x: float((x[x > 1.0] - 1.0).mean()) if (x > 1.0).sum() > 50 else np.nan
idx = np.floor((pos - EXT0) / 4.0).astype(int)
ok = np.all((idx >= 0) & (idx < np.array(TL.shape)), 1)
lab = np.full(len(d), "Out", dtype=object); lab[ok] = TL[idx[ok, 0], idx[ok, 1], idx[ok, 2]]
out = {"E0": {}, "E1": {}, "E2": {}}
for t, c in CFG.items():
    out["E0"][t] = est(d[mat == t]); out["E1"][t] = est(d[lab == t])
# E2: 8 mm voxels, modal label
n8 = (70, 35, 10); i8 = np.floor((pos - EXT0) / 8.0).astype(int)
ok8 = np.all((i8 >= 0) & (i8 < np.array(n8)), 1) & (d > 1.0)
fl = np.ravel_multi_index(i8[ok8].T, n8)
s = np.bincount(fl, weights=d[ok8] - 1.0, minlength=np.prod(n8)).reshape(n8); c = np.bincount(fl, minlength=np.prod(n8)).reshape(n8)
modal = np.empty(n8, dtype=object)
for i in range(n8[0]):
    for j in range(n8[1]):
        for k in range(n8[2]):
            b = TL[2*i:2*i+2, 2*j:2*j+2, 2*k:2*k+2]; v, cc = np.unique(b, return_counts=True); modal[i, j, k] = v[cc.argmax()]
valid = c >= 30
for t in CFG:
    sel = (modal == t) & valid
    out["E2"][t] = float((s[sel]).sum() / c[sel].sum()) if sel.sum() else np.nan
print(f"{'tissue':10s} {'cfg':>5} {'E0 material':>12} {'E1 atlas4mm':>12} {'E2 8mm vox':>12}")
for t, cfg in CFG.items():
    print(f"{t:10s} {cfg:5.2f} " + " ".join(f"{out[k][t]:7.3f}({100*(out[k][t]-cfg)/cfg:+5.1f}%)" for k in ("E0", "E1", "E2")))
sb = lab == "SpineBone"
comp = Counter(mat[sb]); tail = Counter(mat[sb & (d > 1.0)])
print("\natlas-SpineBone cells: annihilation material shares (all / tau>1 ns tail)")
for m, n in comp.most_common(8):
    print(f"  {m:10s} {100*n/sb.sum():5.1f}%   {100*tail.get(m,0)/max((sb&(d>1)).sum(),1):5.1f}%")
out["spinebone_composition"] = {m: [n / int(sb.sum()), tail.get(m, 0) / int((sb & (d > 1)).sum())] for m, n in comp.most_common(8)}
json.dump(out, open(RES / "2026-09-30_referee2/vipman_oracle_mix.json", "w"), indent=1, default=float)
