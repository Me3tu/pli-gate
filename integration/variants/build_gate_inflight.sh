#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# v1.1.0: locations from the environment; preparing SRC is described in integration/PATCH_NOTES.md
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${PLI_ENV_SCRIPT:-$PKG/simulation/runner/env_example.sh}" >/dev/null 2>&1
SRC=${GATE_INFLIGHT_SRC:?set GATE_INFLIGHT_SRC to a GATE 9.4.1 tree prepared with install_into_gate.sh and gate-inflight-stock.patch}
BLD=${GATE_INFLIGHT_BUILD:-$PWD/build-gate-inflight}
INS=${GATE_INFLIGHT_INSTALL:-$PWD/install-gate-inflight}
mkdir -p $BLD && cd $BLD
date +%T > build.start
cmake -S $SRC -B $BLD -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INS \
  -DCMAKE_CXX_COMPILER=/usr/bin/c++ -DCMAKE_C_COMPILER=/usr/bin/cc \
  -DGeant4_DIR=$GEANT4_INSTALL/lib/cmake/Geant4 -DROOT_DIR=$ROOT_INSTALL/cmake \
  -DGATE_USE_OPTICAL=ON -DGATE_USE_RTK=OFF -DGATE_USE_ITK=OFF -DGATE_USE_LMF=OFF -DGATE_USE_DAVIS=OFF \
  -DGATE_USE_TORCH=OFF -DGATE_USE_XRAYLIB=OFF -DGATE_USE_SYSTEM_CLHEP=OFF -DGATE_COMPILE_GATEDIGIT=OFF \
  -DGATE_COMPILE_WITH_NATIVE=OFF > cmake.log 2>&1 || { echo CMAKE_FAIL > build.status; exit 1; }
make -j28 > make.log 2>&1 || { echo MAKE_FAIL > build.status; exit 1; }
make install > install.log 2>&1 || { echo INSTALL_FAIL > build.status; exit 1; }
date +%T > build.end
echo BUILD_OK > build.status
