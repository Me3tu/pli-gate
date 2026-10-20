# How the PLI physics enters GATE

Two mechanisms, one trivial and one small. A third section describes the
in-flight annihilation variant used for one sensitivity arm of the paper.

## 1. Source files: a plain copy, no build-system change

GATE's top-level `CMakeLists.txt` gathers sources by glob:

```cmake
FILE(GLOB sources
  ...
  ${PROJECT_SOURCE_DIR}/source/physics/src/*.cc
  ...)
FILE(GLOB headers
  ...
  ${PROJECT_SOURCE_DIR}/source/physics/include/*.hh
  ...)
```

So five `.cc` files copied into `source/physics/src/` and five `.hh` files
copied into `source/physics/include/` are compiled into GATE automatically.
There is nothing to add to any CMake file.

One consequence worth knowing: **any** file ending in `.cc` in that directory
is compiled. Do not leave editor backups named `something.cc` there. Backups
ending in `.backup` or `.bak` are ignored because the glob matches `*.cc`.

`PLISteppingAction` is *not* installed into GATE. It belongs to the optional
standalone library and is not part of the GATE build.

## 2. The patch: one hook plus portability fixes

`gate-9.4.1-pli.patch` applies to GATE 9.4.1 at upstream commit
`2d03a89bfec37c21ec69d196b4876d5e3bf5ed83`. It touches 11 files. GATE was
patched, Geant4 was not. The patch is byte-identical to the difference
between that commit and the GATE tree that built the production binary.

**The one functional change** is in `source/general/src/GateRunManager.cc`: it
lets a macro select the physics list by name,

```
/gate/physics/addPhysicsList PLI_Physics
```

by constructing `PLIPhysicsList` directly instead of asking
`G4PhysListFactory` for a name it does not know. Without this hook the physics
compiles but cannot be selected.

**The other ten changes are portability and robustness**, not physics, and are
needed to build and run this GATE version on a current toolchain:

| files | change | why |
|---|---|---|
| `source/io/**` (6 files) | add `#include <cstdint>` | GCC 15 no longer provides these types transitively |
| `source/externals/itk-mhd/CMakeLists.txt` | `cmake_minimum_required` 2.4 → 3.5 | CMake 4 refuses compatibility below 3.5 |
| `GateOutputMgr.cc`, `GateToRoot.cc` | null / empty-container checks | avoid a crash when an output module is absent or no tree was filled |
| `GateOutputMgr.cc`, `GateUserActions.cc` | disable four `/vis/...` calls | these crash under Geant4 11.3.2 when no visualisation is active |

None of the ten alters simulated physics. They are separated here so that a
reader can see the physics change is exactly one hook.

## 3. Variant: stock in-flight annihilation kept (`variants/`)

One sensitivity arm of the paper (the in-flight arm, binary cc552a68) used a
physics list that keeps Geant4's own `annihil` process for positrons that
annihilate in flight and deactivates only its at-rest part, so that every
positron that stops is still handled by the PLI processes. The switch is the
environment variable `PLI_INFLIGHT_STOCK=1`. Unset, or with any other value,
the variant behaves exactly like the shipped list.

| file | what it is |
|---|---|
| `variants/gate-inflight-stock.patch` | the whole source difference, one hunk against `custom_physics/src/PLIPhysicsList.cc`, made with `diff -u` |
| `variants/build_gate_inflight.sh` | the build script of the variant binary, with the production configuration |
| `variants/pli_inflight_test/` | the physics-level test: a standalone Geant4 program that compares PLI, variant and stock annihilation for 44Sc or 124I positrons (`CMakeLists.txt`, `pli_inflight_test.cc`, `run_physics_level.sh`) |

Building the variant GATE:

```bash
cp -r $GATE_SRC $GATE_INFLIGHT_SRC                       # a second, clean GATE 9.4.1 tree
$PKG/integration/install_into_gate.sh $GATE_INFLIGHT_SRC
patch $GATE_INFLIGHT_SRC/source/physics/src/PLIPhysicsList.cc \
      < $PKG/integration/variants/gate-inflight-stock.patch
GATE_INFLIGHT_SRC=$GATE_INFLIGHT_SRC $PKG/integration/variants/build_gate_inflight.sh
```

The build script reads `GEANT4_INSTALL` and `ROOT_INSTALL` from
`simulation/runner/env_example.sh` (or from `PLI_ENV_SCRIPT`), and installs to
`GATE_INFLIGHT_INSTALL` (default `./install-gate-inflight`). Point
`GATE_INFLIGHT_BIN_DIR` at its `bin/` before `simulation/runs/inflight_arm.sh`.

The physics-level test applies the patch itself, in its build tree, at
configure time:

```bash
cmake -S $PKG/integration/variants/pli_inflight_test -B build -DGeant4_DIR=$GEANT4_INSTALL/lib/cmake/Geant4
cmake --build build
```

Checked for this release: the shipped `PLIPhysicsList.cc` with this patch
compiles, with the production flags, to an object file byte-identical to the
one in the cc552a68 build. The test built from this folder reproduces the
output of the production test program byte for byte (3,000 44Sc positrons,
seed 401, in each of the modes `pli`, `hybrid` and `stock`).

## Reversing

The patch is a normal unified diff: `git apply -R gate-9.4.1-pli.patch`, or
`patch -p1 -R < gate-9.4.1-pli.patch`. Removing the ten PLI files from
`source/physics/` returns GATE to stock.
