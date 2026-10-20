#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
"""R2c materials: database v8 with I3 = 0.2625 for every o-Ps-intensity entry (P_oPs = 0.2625 in every material,
matching the total-formation structural entries at I3 = 0.35). Writes r2c_equal_yield_materials.json here, and
the released database is read only."""
import json, hashlib
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]          # v1.1.0: release root
HERE = Path(__file__).resolve().parent
src = REPO / 'materials/pli_materials.json'
d = json.loads(src.read_text())
n = 0
for e in d['materials']:
    if e.get('I3_convention') == 'ops_intensity':
        e['I3_r2c_original'] = e['I3']; e['I3'] = 0.2625; n += 1
    else:
        assert e.get('I3_convention') == 'total_formation' and abs(e['I3'] - 0.35) < 1e-12, e['name']
d.setdefault('_r2c', {})['note'] = ('R2c matched-yield sensitivity (2026-10-08): I3 = 0.2625 for every o-Ps-intensity '
                                   'entry; derived from %s sha256 %s' % (src.name, hashlib.sha256(src.read_bytes()).hexdigest()))
out = HERE / 'r2c_equal_yield_materials.json'
out.write_text(json.dumps(d, indent=1))
print('entries changed', n, '->', out, hashlib.sha256(out.read_bytes()).hexdigest()[:12])
