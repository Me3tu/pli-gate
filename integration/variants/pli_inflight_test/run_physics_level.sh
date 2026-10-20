#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# v1.1.0: expects the CMake build of this test in ./build (cmake -S . -B build), and
# output goes to $PLI_OUT/feedback3_inflight/physics.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${GEANT4_INSTALL:?set GEANT4_INSTALL}/bin/geant4.sh"
export PLI_MATERIALS_JSON=${PLI_MATERIALS_JSON:-$HERE/../../../materials/pli_materials.json}
PHYS=${PLI_OUT:?set PLI_OUT}/feedback3_inflight/physics; mkdir -p "$PHYS"
cd "$HERE"
for iso in sc44 i124; do
  for m in pli hybrid stock; do
    seed=$([ $iso = sc44 ] && echo 401 || echo 501)
    ( nice -n 10 ./build/pli_inflight_test $m $iso 1000000 1.0 $seed $PHYS/${iso}_${m}.csv > $PHYS/${iso}_${m}.log 2>&1; echo "$(date +%T) $iso $m rc=$?" >> $PHYS/physics_level.log ) &
  done
done
wait
echo "$(date +%T) ALL_DONE" >> $PHYS/physics_level.log
