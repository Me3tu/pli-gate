#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-only
# Copyright (C) 2026 Dhritman Dhar and German Cancer Research Center (DKFZ)
#
# Install the PLI physics into a GATE 9.4.1 source tree.
#
# GATE's top-level CMakeLists.txt collects its sources with
#
#     FILE(GLOB sources ... ${PROJECT_SOURCE_DIR}/source/physics/src/*.cc ...)
#     FILE(GLOB headers ... ${PROJECT_SOURCE_DIR}/source/physics/include/*.hh ...)
#
# so a source file placed in those two directories is compiled into GATE with
# no build-system change at all. That is the whole integration: this script
# copies files, and the patch adds the one hook that lets a macro select the
# physics list by name.
#
# Usage:  ./install_into_gate.sh /path/to/gate-9.4.1
#
set -euo pipefail

GATE_SRC="${1:-}"
if [[ -z "$GATE_SRC" ]]; then
    echo "usage: $0 /path/to/gate-9.4.1" >&2; exit 2
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/.." && pwd)"

if [[ ! -d "$GATE_SRC/source/physics/src" ]]; then
    echo "ERROR: $GATE_SRC does not look like a GATE source tree" >&2
    echo "       (expected $GATE_SRC/source/physics/src)" >&2; exit 1
fi

SRCS=(G4PositroniumAnnihilation G4PsAtRestProcess PLIPhysicsList
      PLIMaterialLoader PLIOrePowell)
HDRS=(G4PositroniumAnnihilation G4PsAtRestProcess PLIPhysicsList
      PLIMaterialLoader PLIOrePowell)

echo "==> installing PLI physics into $GATE_SRC"
for f in "${SRCS[@]}"; do
    cp -v "$PKG/custom_physics/src/$f.cc" "$GATE_SRC/source/physics/src/$f.cc"
done
for f in "${HDRS[@]}"; do
    cp -v "$PKG/custom_physics/include/$f.hh" "$GATE_SRC/source/physics/include/$f.hh"
done

# nlohmann/json is required by PLIMaterialLoader. GATE's include path already
# covers source/physics/include, so vendoring the single header there keeps the
# build self-contained and needs no system package.
NL_DST="$GATE_SRC/source/physics/include/nlohmann"
if [[ -f "$PKG/custom_physics/third_party/nlohmann/json.hpp" ]]; then
    mkdir -p "$NL_DST"
    cp -v "$PKG/custom_physics/third_party/nlohmann/json.hpp" "$NL_DST/json.hpp"
else
    echo "    nlohmann/json.hpp not bundled; relying on a system installation"
fi

echo "==> applying the GATE patch"
if git -C "$GATE_SRC" rev-parse --git-dir >/dev/null 2>&1; then
    git -C "$GATE_SRC" apply --check "$HERE/gate-9.4.1-pli.patch" 2>/dev/null \
        && git -C "$GATE_SRC" apply "$HERE/gate-9.4.1-pli.patch" \
        && echo "    patch applied" \
        || echo "    patch already applied (or does not apply cleanly) - see PATCH_NOTES.md"
else
    patch -p1 -d "$GATE_SRC" -N --dry-run < "$HERE/gate-9.4.1-pli.patch" >/dev/null 2>&1 \
        && patch -p1 -d "$GATE_SRC" -N < "$HERE/gate-9.4.1-pli.patch" \
        && echo "    patch applied" \
        || echo "    patch already applied (or does not apply cleanly) - see PATCH_NOTES.md"
fi

echo
echo "==> done. Next: configure and build GATE against your Geant4 and ROOT,"
echo "    then run examples/quickstart/. See the README."
