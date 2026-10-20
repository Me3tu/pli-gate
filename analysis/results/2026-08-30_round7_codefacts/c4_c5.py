#!/usr/bin/env python3
"""C4 tie-rate + C5 t_hat_ann - t_true, both on the real truth-tagged stream.

C4: how often do two in-511 singles share an EXACT timestamp, so that the
    unstable sort in M3_clinical_extractor_randoms.py:164 could reorder them?
C5: bias and spread of the TOF-reconstructed annihilation TIME.
    truth uses the UNSMEARED time_ns and the truth source position:
        t_true = mean_i [ t_true_i - |g_i - src|/c ]
    estimate uses the SMEARED time and the TOF annihilation point:
        t_hat  = mean_i [ t_smear_i - |g_i - p_ann|/c ]
    Averaging over the two photons is deliberate: src is the DECAY vertex, not
    the annihilation vertex, so the per-photon path error from positron range
    is +d for one photon and -d for the other and cancels in the mean.
"""
import os, sys, json
import numpy as np, h5py
from pathlib import Path
SRC=str(Path(os.environ.get("PLI_RAW", "<missing>")) / "2026-06-10_brain_tumor_1B_truetime/sim1_singles_truetime.h5")
if not Path(SRC).exists(): raise SystemExit(f"{SRC}: raw data not included in the release, see analysis/README.md")
C=299.792458  # mm/ns
W_PAIR, SIG_T = 2.0, 0.0909
N=8_000_000
with h5py.File(SRC,"r") as f:
    g=f["singles"]
    sh=g["shard_id"][:N]
    end=int(np.flatnonzero(sh!=sh[-1]).max())+1
    t=g["time_ns"][:end]; i5=g["in_511"][:end].astype(bool)
    gx,gy,gz=g["gx"][:end],g["gy"][:end],g["gz"][:end]
    sx,sy,sz=g["sx"][:end],g["sy"][:end],g["sz"][:end]
    shard=sh[:end]
print(f"  slice {end:,} singles, {len(np.unique(shard))} shard(s), in_511 {int(i5.sum()):,}")

rng=np.random.default_rng(20260830)
t_sm=t+rng.normal(0.0,SIG_T,size=len(t))

res={}
for lab,tt in (("TRUE (unsmeared)",t),("SMEARED",t_sm)):
    idx=np.flatnonzero(i5); ts=tt[idx]
    o=np.argsort(ts,kind="stable"); idx=idx[o]; ts=ts[o]
    d=np.diff(ts)
    res[lab]=dict(n=len(ts),ties=int((d==0).sum()),pairs=int((d<=W_PAIR).sum()))
    print(f"  {lab:<18} n={len(ts):,}  exact ties={int((d==0).sum()):,}  "
          f"({100*(d==0).mean():.4f} %)  pairs(dt<=2ns)={int((d<=W_PAIR).sum()):,}")

# ---- C4: can one single be in two pairs? measure it -----------------------
idx=np.flatnonzero(i5); ts=t_sm[idx]
o=np.argsort(ts,kind="stable"); idx=idx[o]; ts=ts[o]
ok=np.diff(ts)<=W_PAIR
pa,pb=idx[:-1][ok],idx[1:][ok]
shared=int(np.intersect1d(pa,pb).size)
both=np.concatenate([pa,pb]); u,c=np.unique(both,return_counts=True)
print(f"\n  C4  pairs formed              {len(pa):,}")
print(f"      singles used in >1 pair   {int((c>1).sum()):,} "
      f"({100*(c>1).mean():.3f} % of participating singles)")
print(f"      max pairs one single joins {int(c.max())}")
print(f"      (a single that is pb of one pair and pa of the next: {shared:,})")

# ---- C5: TOF annihilation-time error --------------------------------------
p1=np.c_[gx[pa],gy[pa],gz[pa]]; p2=np.c_[gx[pb],gy[pb],gz[pb]]
src=np.c_[sx[pa],sy[pa],sz[pa]]
mp=0.5*(p1+p2); lor=p2-p1
lu=lor/np.maximum(np.linalg.norm(lor,axis=1,keepdims=True),1e-9)
pann=mp+0.5*C*(t_sm[pa]-t_sm[pb])[:,None]*lu
d1h=np.linalg.norm(p1-pann,axis=1); d2h=np.linalg.norm(p2-pann,axis=1)
d1t=np.linalg.norm(p1-src,axis=1);  d2t=np.linalg.norm(p2-src,axis=1)
t_hat =0.5*((t_sm[pa]-d1h/C)+(t_sm[pb]-d2h/C))
t_true=0.5*((t[pa]  -d1t/C)+(t[pb]  -d2t/C))
err=(t_hat-t_true)*1000.0   # ps
# restrict to TRUE coincidences (same decay): both singles share the source pos
same=np.all(np.isclose(src,np.c_[sx[pb],sy[pb],sz[pb]],atol=1e-6),axis=1)
print(f"\n  C5  pairs {len(err):,}; same-origin (true coincidences) {int(same.sum()):,} "
      f"({100*same.mean():.1f} %)")
for lab,m in (("all pairs",np.ones(len(err),bool)),("true coincidences only",same)):
    e=err[m]; e=e[np.isfinite(e)]
    print(f"      {lab:<24} mean {np.mean(e):+8.2f} ps   sd {np.std(e,ddof=1):8.2f} ps   "
          f"median {np.median(e):+8.2f} ps   IQR {np.subtract(*np.percentile(e,[75,25])):.2f} ps")
e=err[same]; e=e[np.isfinite(e)]
res["c5"]=dict(mean_ps=float(np.mean(e)),sd_ps=float(np.std(e,ddof=1)),
               median_ps=float(np.median(e)),n=int(len(e)))
(Path(os.environ.get("PLI_RES") or Path(__file__).resolve().parents[1]) / "2026-08-30_round7_codefacts/c4_c5.json").write_text(json.dumps(res,indent=2))
