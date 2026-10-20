# Building the dependencies from source

The README assumes you already have GATE 9.4.1, Geant4 11.3.2 and ROOT 6.36.04.
If you do not, this is the path from nothing to a working install, with the
exact options used to validate this package. Everything lands under one
prefix and nothing is installed system-wide.

```bash
export SC=$HOME/pli-stack          # one directory holds everything
mkdir -p $SC/{src,build,install,logs}
```

Get the three sources into `$SC/src/{geant4,root,gate}`. GATE must be at
upstream commit `2d03a89bfec37c21ec69d196b4876d5e3bf5ed83`; see
`integration/gate-base-commit.txt`.

## 1. Geant4 11.3.2

```bash
cmake -B $SC/build/geant4 -S $SC/src/geant4 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$SC/install/geant4 \
  -DGEANT4_INSTALL_DATA=ON \
  -DGEANT4_BUILD_MULTITHREADED=OFF \
  -DGEANT4_USE_GDML=ON \
  -DGEANT4_USE_QT=OFF \
  -DGEANT4_USE_OPENGL_X11=OFF \
  -DGEANT4_USE_RAYTRACER_X11=OFF
cmake --build $SC/build/geant4 -j $(nproc)
cmake --install $SC/build/geant4
```

`GEANT4_USE_GDML=ON` needs xerces-c. `GEANT4_INSTALL_DATA=ON` downloads about
2.1 GB of physics datasets; if you already have them, set it `OFF` and point
`GEANT4_INSTALL_DATADIR` at the existing directory instead.

**Single-threaded on purpose.** `GEANT4_BUILD_MULTITHREADED=OFF` matches the
configuration this package was validated in.

## 2. ROOT 6.36.04

GATE calls `FIND_PACKAGE(ROOT REQUIRED)` with no components and links
`${ROOT_LIBRARIES}`, so it needs Core, RIO, Tree, Hist, MathCore and Matrix.
Everything else can be switched off, which cuts the build substantially.

```bash
cmake -B $SC/build/root -S $SC/src/root \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$SC/install/root \
  -Dxrootd=OFF -Dmysql=OFF -Dpgsql=OFF -Doracle=OFF -Dodbc=OFF \
  -Dpythia6=OFF -Dpythia8=OFF -Dfftw3=OFF -Droofit=OFF -Dtmva=OFF \
  -Dhttp=OFF -Ddavix=OFF -Dwebgui=OFF -Dgviz=OFF -Dopengl=OFF \
  -Dx11=OFF -Dasimage=OFF -Dpyroot=OFF -Dtpython=OFF \
  -Dspectrum=OFF -Dunuran=OFF -Dvdt=OFF -Dxml=ON \
  -Dtesting=OFF -Droottest=OFF
cmake --build $SC/build/root -j $(nproc)
cmake --install $SC/build/root
```

**`-Dtpython=OFF` is not optional if you set `-Dpyroot=OFF`.** ROOT enables
`tpython` by default and its configure aborts with

```
CMake Error at cmake/modules/RootBuildOptions.cmake:456 (message):
  -Dtpython=ON requires -Dpyroot=ON)
```

Turning one off without the other is the first thing that goes wrong here.

## 3. The environment

Source this before building GATE and before every run. It is the whole
environment — nothing else is needed and nothing system-wide is used.

```bash
source $SC/install/geant4/bin/geant4.sh
export ROOTSYS=$SC/install/root
export PATH=$SC/install/gate/bin:$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$SC/install/geant4/lib:${LD_LIBRARY_PATH:-}
export PLI_MATERIALS_JSON=<path-to-this-package>/materials/pli_materials.json
```

The run scripts of `simulation/` source `simulation/runner/env_example.sh`,
which expects this layout. Set `PLI_STACK=$SC/install` and it finds Geant4,
ROOT and GATE under `$PLI_STACK/{geant4,root,gate}`.

## 4. GATE and this package

Now follow the README from step 1, pointing at your own prefixes:

```bash
$PKG/integration/install_into_gate.sh $SC/src/gate
cmake -B $SC/build/gate -S $SC/src/gate \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$SC/install/gate \
      -DGeant4_DIR=$SC/install/geant4/lib/cmake/Geant4 \
      -DROOT_DIR=$ROOTSYS/cmake
cmake --build $SC/build/gate -j $(nproc)
cmake --install $SC/build/gate
```

Then the quickstart and `tests/run_tests.sh` as the README describes.

## Measured

Timings from a validated from-source bootstrap on a 32-core desktop, Arch
Linux, GCC 15.2.1, CMake 4.2.1.

| stage | wall | cores | notes |
|---|---|---|---|
| Geant4 configure + build + install | 275 s | 16 | 33 libraries, GDML on, single-threaded |
| ROOT configure + build + install | 906 s | 16 | reduced feature set as above |
| GATE + this package | 90 s | 24 | includes quickstart and the test suite |
| **total, nothing to working install** | **~21 min** | | |

Disk:

| | |
|---|---|
| Geant4 install (2.1 GB of that is physics datasets) | 2.3 GB |
| ROOT install | 433 MB |
| GATE install | 13 MB |
| build trees (removable afterwards) | 2.0 GB |
| **total** | **5.4 GB** |

## What this configuration was validated to do

Built exactly as above and then checked, with the whole toolchain freshly
compiled and nothing system-wide on the path:

* quickstart recovers **tau_oPs = 1.9072 +/- 0.0107 ns** against a configured
  1.9215 ns
* formation channels match theory: o-Ps share 0.2694 against I3 = 0.27, p-Ps
  0.0904 against I3/3, direct 0.6402 against 0.6400, three-photon share of o-Ps
  0.01345 against tau/142.05 = 0.01353, zero envelope violations
* all nine tests pass
* an undeclared material aborts the run and names both the material and the
  escape hatch

The run produced **200,703 annihilations and 31,894 above the censor**, the
same counts as the reference run in `examples/quickstart/expected_output.json`.
Rebuilding the same Linux stack reproduces those counts exactly.

Across platforms the counts drift and the lifetime does not. A macOS/clang
build of the same GATE, Geant4 and package gave 200,683 annihilations
(0.010 % fewer) and 31,880 above the censor (0.044 % fewer), yet returned
tau_oPs = 1.9071 ns against 1.9072 ns here, agreeing to 0.1 ps. Geant4 does not
guarantee bit-identical sampling across compilers and standard libraries, so
treat the counts as toolchain-dependent and the lifetime as the quantity that
travels. `expected_output.json` states the same thing in its own note.
