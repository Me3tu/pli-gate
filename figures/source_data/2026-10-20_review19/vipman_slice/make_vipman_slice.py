"""Input of Fig. S7 panel (a): the central two 4 mm slices of the corrected-yield VIPman MLEM
reconstruction, in float32 (2026-10-20, review round 19, for code release v1.1.0).

Replaces, for the release, the 8-bit transfer 2026-10-19_review16/vipman_recon_O4/slice_q8.b64 from
which the panel of the submitted figure was drawn (step 7.24 a.u., values above 1846.2 clipped).

Source (../../2026-10-19_review16/vipman_recon_O4/, simple_recon single, 4 mm voxels, 10 MLEM
iterations, grid 140 x 70 x 20, the settings of recon_phase3q):
  activity_x.npy    float32 (140, 70, 20)
  support_mask.npy  bool    (140, 70, 20)

Output (next to this script): vipman_O4_slice.npz, written with numpy.savez_compressed
  act_z  float32 (140, 70)     mean of z slices 10 and 11 of activity_x.npy, no mask, no rescaling
  sup2   bool    (140, 70, 2)  support_mask.npy[:, :, [10, 11]]
f05_vipman.py draws act_z only where sup2.any(axis=2) and asserts that sup2 equals the support mask
of recon_phase3q in the same slices.

Run: python3 make_vipman_slice.py   (prints input hashes and hashes of the array bytes. The npz
bytes themselves carry zip timestamps, so compare the printed array hashes across reruns).
"""
import hashlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC = HERE.parent.parent / '2026-10-19_review16' / 'vipman_recon_O4'
Z4 = [10, 11]
OUT = HERE / 'vipman_O4_slice.npz'


def sha_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sha_arr(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


act = np.load(SRC / 'activity_x.npy')
sup = np.load(SRC / 'support_mask.npy')
assert act.shape == (140, 70, 20) and act.dtype == np.float32, (act.shape, act.dtype)
assert sup.shape == (140, 70, 20) and sup.dtype == np.bool_, (sup.shape, sup.dtype)

act_z = act[:, :, Z4].mean(axis=2).astype(np.float32)   # same operation as f05_vipman.py on recon_phase3q
sup2 = np.ascontiguousarray(sup[:, :, Z4])
np.savez_compressed(OUT, act_z=act_z, sup2=sup2)

m = sup2.any(axis=2)
shown = np.where(m, act_z, np.nan).T
print('source activity_x.npy  ', sha_file(SRC / 'activity_x.npy'))
print('source support_mask.npy', sha_file(SRC / 'support_mask.npy'))
print('act_z', act_z.dtype, act_z.shape, 'sha256(bytes)', sha_arr(act_z))
print('sup2 ', sup2.dtype, sup2.shape, 'sha256(bytes)', sha_arr(sup2))
print('support voxels: sup2', int(sup2.sum()), 'sup2.any(axis=2)', int(m.sum()))
print('inside support: min %.4f max %.4f 99th percentile (colour-scale top of Fig. S7a) %.4f'
      % (np.nanmin(shown), np.nanmax(shown), np.nanpercentile(shown, 99)))
print('wrote', OUT, OUT.stat().st_size, 'bytes, sha256', sha_file(OUT))
