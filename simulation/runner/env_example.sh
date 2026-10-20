#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
# Single-threaded environment for PLI validation

# v1.1.0: install locations come from the environment; the defaults follow the
# one-prefix layout of docs/BOOTSTRAP.md ($PLI_STACK/{geant4,root,gate}).
export PLI_STACK=${PLI_STACK:-${HOME}/pli-stack/install}

# Single-threaded installations
export GEANT4_INSTALL=${GEANT4_INSTALL:-${PLI_STACK}/geant4}
export ROOT_INSTALL=${ROOT_INSTALL:-${PLI_STACK}/root}
export GATE_INSTALL=${GATE_INSTALL:-${PLI_STACK}/gate}

# Source Geant4 environment
source ${GEANT4_INSTALL}/bin/geant4.sh

# Set up ROOT manually (thisroot.sh is deprecated)
export ROOTSYS=${ROOT_INSTALL}
export PATH=${ROOT_INSTALL}/bin:${PATH}
export LD_LIBRARY_PATH=${ROOT_INSTALL}/lib:${LD_LIBRARY_PATH:-}
export PYTHONPATH=${ROOT_INSTALL}/lib:${PYTHONPATH:-}

# Add GATE to PATH (prepend to override system GATE if any)
export PATH=${GATE_INSTALL}/bin:${PATH}

# Ensure LD_LIBRARY_PATH includes Geant4
export LD_LIBRARY_PATH=${GEANT4_INSTALL}/lib:${LD_LIBRARY_PATH}

echo "✓ Geant4 environment loaded (ST)"
echo "✓ ROOT environment loaded"
echo "✓ GATE added to PATH"
echo ""
echo "Software stack ready (single-threaded):"
echo "  Geant4: $(geant4-config --version)"
echo "  ROOT:   $(root-config --version)"
echo "  GATE:   Gate version is 9.4.1"
