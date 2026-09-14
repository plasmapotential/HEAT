# AGENTS.md — HEAT

Guidance for AI agents (and humans) working in this repository. Part 1 covers the
repo itself: running, testing, releasing, and how the code is organized. Part 2 is the
runbook for preparing and executing a HEAT simulation case.

## What is HEAT

The **Heat flux Engineering Analysis Toolkit (HEAT)** is a Python suite for predicting heat flux incident on plasma-facing components (PFCs) in tokamaks. It combines CAD geometry, MHD equilibria, and multiple heat flux models (optical, ion gyro orbit, photon radiation, filaments, runaway electrons, 3D fields) into one framework. Developed by Tom Looby at Commonwealth Fusion Systems; used to design SPARC PFCs.

---

# Part 1 — Working in the repository

## Running HEAT

HEAT runs inside Docker. The published image is `plasmapotential/heat:<tag>`; the current tag is `HEAT_IMAGE_TAG` in `.github/workflows/integration-tests.yml` and the `image:` line of `docker/docker-compose.yml`. Substitute it for `<tag>` in every command below.

**Start the GUI (web app on localhost:8050):**
```bash
cd docker && docker compose up
```

**TUI/batch mode (inside container or from compose):**
```bash
docker run --rm -v "$(pwd):/root/source/HEAT" plasmapotential/heat:<tag> \
  --m t --f /root/source/HEAT/tests/integrationTests/nstxuTestCase/batchFile_optical.dat
```

**Interactive shell in container:**
```bash
docker compose run --entrypoint "" HEAT /bin/bash
```

The `docker/docker-compose.yml` already mounts `~/HEAT` for data and the local repo source into the container, so local code changes are picked up immediately without rebuilding the image.

## Running Tests

All tests run inside the Docker container against the published image (they rely on the full compiled environment: FreeCAD, OpenFOAM, Open3D, Mitsuba).

**Smoke test (sanity check the mount):**
```bash
docker run --rm -v "$(pwd):/root/source/HEAT" --entrypoint "" \
  plasmapotential/heat:<tag> \
  python3 /root/source/HEAT/tests/integrationTests/ciTest.py
```

**Single integration test case (e.g. optical):**
```bash
docker run --rm -v "$(pwd):/root/source/HEAT" \
  plasmapotential/heat:<tag> \
  --m t --f /root/source/HEAT/tests/integrationTests/nstxuTestCase/batchFile_optical.dat
```

Available batch files in `tests/integrationTests/nstxuTestCase/`:
- `batchFile_optical.dat` — optical heat flux
- `batchFile_optical_elmer.dat` — optical + Elmer FEM thermal solve
- `batchFile_gyro.dat` — ion gyro orbit
- `batchFile_rad.dat` — photon radiation (also has golden assertions)
- `batchFile_rzq.dat` — R,Z,q|| profile from CSV
- `batchFile_optical_BYOM.dat` — bring your own mesh (STL input)

**Photon radiation golden checks:**
```bash
python3 tests/integrationTests/verify_nstxu_hf_rad_goldens.py \
  --workspace "$(pwd)" --docker-image plasmapotential/heat:<tag>

# or via pytest:
pytest tests/integrationTests/test_nstxu_hf_rad_goldens.py -v
```

CI runs all of the above automatically on push/PR to `main` (`.github/workflows/integration-tests.yml`).

**Updating goldens** when physics intentionally change: re-run the rad verifier, copy the printed "Parsed metrics" into `tests/integrationTests/nstxuTestCase/nstxu_hf_rad_goldens.json`.

## Release

```bash
./scripts/release.sh <IMAGE_TAG> [HEAT_REF] [--build]
# e.g.:
./scripts/release.sh v4.2.8 v4.3 --build
```

This updates `HEAT_IMAGE_TAG` in CI and the docker-compose image tags, optionally builds the image, then prompts you to push and open a PR to `main`.

## Architecture

### Entry point and modes

`source/launchHEAT.py` is the single entry point. It reads the `runMode` environment variable (`docker` or `local`), sets up all paths and `sys.path` entries for external tools (FreeCAD, ParaView, OpenFOAM, EFIT, Open3D), then hands off to either:
- `dashGUI.py` — Plotly Dash web application (GUI mode, `--m g`)
- `terminalUI.py` — batch/terminal mode (`--m t --f batchFile.dat`)

### Engine and physics modules

`engineClass.engineObj` is the central orchestrator. It owns one instance of every physics module and coordinates the time-stepping loop. All modules are instantiated in `engineClass.initializeEveryone()`:

| Engine attribute | Class | Role |
|---|---|---|
| `ENG.MHD` | `MHDClass.MHD` | Reads GEQDSK/gfile equilibria (via EFIT class), field-line mapping |
| `ENG.CAD` | `CADClass.CAD` | Loads STEP or STL geometry via FreeCAD, meshes PFCs |
| `ENG.HF` | `heatfluxClass.heatFlux` | Optical heat flux (Eich profile, multiExp, tophat, qFile) |
| `ENG.GYRO` | `gyroClass.GYRO` | Ion gyro-orbit heat flux tracing |
| `ENG.RAD` | `radClass.RAD` | Photon/radiation heat flux (uses Mitsuba for ray tracing) |
| `ENG.FIL` | `filamentClass.filament` | ELM filament heat and particle fluxes |
| `ENG.RE` | `runawayClass.Runaways` | Runaway electron module |
| `ENG.plasma3D` | `plasma3DClass.plasma3D` | 3D field perturbations via MAFOT/laminar |
| `ENG.OF` | `openFOAMclass.OpenFOAM` | OpenFOAM thermal conduction solver |
| `ENG.FEM` | `elmerClass.FEM` | Elmer FEM thermal solver (via gmsh) |
| `ENG.IO` | `ioClass.IO_HEAT` | Output: VTP meshes, point clouds, CSV |

`rayTracerClass` (extracted from `pfcClass.py`) provides the shared ray–mesh intersection kernels used by RAD, FIL, RE, and PFC shadow detection (wraps Open3D and Mitsuba).

### PFC object — the central data structure

`pfcClass.PFC` inherits from `rayTracerClass.shadowKernels`. One PFC object is created per tile defined in the PFC input CSV. It holds the CAD mesh data, the MHD equilibrium data mapped onto that tile, and the heat flux profile parameters. Shadow detection (which field lines are blocked by other tiles) is done per PFC at each timestep.

### Input pipeline

HEAT is configured entirely via CSV input files. Each module has an `allowed_class_vars()` method listing recognized parameter names; loading an input CSV populates the corresponding module's attributes. For batch/TUI mode, a `batchFile.dat` lists one HEAT job per row with columns: `MachFlag, Tag, Shot, TimeStep, EQ, CAD, PFC, Input, Output`.

### External dependencies (not pip-installable alone)

- **EFIT class** (ORNL): loaded from `~/source/EFIT/` — reads GEQDSK format equilibria
- **FreeCAD**: for STEP→mesh; path set in `launchHEAT.py`
- **MAFOT**: external binaries (`heatstructure`, `heatlaminar_mpi`) for field-line integration, including 3D perturbed fields. Built inside the image by `docker/buildMAFOT` from MAFOT's `install/make.inc.HEAT`. With `mafot_gpu, True` in the input file HEAT passes `-g` and MAFOT traces on the GPU; the binary only runs on GPU architectures listed in that file's `NVCCFLAGS` (`-gencode` per `sm_XY`), so add an entry there before deploying to a new GPU type. `MHDClass.runMAFOT` raises if a MAFOT call exits non-zero.
- **OpenFOAM**: thermal conduction solver launched as subprocess
- **Mitsuba / drjit**: photon ray tracing in `rayTracerClass.py` and `radClass.py`

All of these are pre-installed in the Docker image. Local development outside Docker requires manually building/installing each one.

### Output

HEAT writes results to `~/HEAT/data/<machine>/<shot>/<timestep>/` (controlled by `heatdata` env var). Each PFC gets a subdirectory. Output formats are controlled by `ioClass` variables (`vtpMeshOut`, `vtpPCOut`, `csvOut`). VTP files are loaded in ParaView for visualization.

### GUIscripts

`source/GUIscripts/` contains plotting helpers used by both GUI and TUI:
- `plotlyGUIplots.py` — Plotly figures served in the Dash GUI
- `plotly2DEQ.py` — equilibrium cross-section plots
- `vtkOpsClass.py` — VTK/VTP file construction
- `meshOpsClass.py` — GLTF/USD mesh export
- `postProcessFunctions.py` — post-processing utilities

### Supported machines

`MachFlag` selects the tokamak: `sparc`, `arc`, `cmod`, `d3d`, `nstx`, `st40`, `step`, `west`, `kstar`, `aug`, `tcv`, `other`. Machine selection in `engineClass.setInitialFiles()` sets CAD and mesh directory paths under `dataPath`.

---

# Part 2 — Preparing and running a HEAT case

## Mental model

- A **HEAT case** is a self-contained directory anywhere on the host, conventionally
  `~/HEATruns/<MACHINE>/<caseName>/`. It contains a `batchFile.dat` plus one
  subdirectory named after the machine flag (e.g. `sparc/`) holding every input file.
- HEAT executes **inside the Docker container** (`plasmapotential/heat:<tag>`). The case
  directory is bind-mounted at `/root/terminal`, so every path in the batch file is
  resolved as `/root/terminal/<MachFlag>/<file>`.
- One batch file row = one (tag, timestep) combination. Rows sharing a **Tag** form one
  simulation; the CAD/PFC/Input/Output columns are only read from the **first row** of
  each tag.
- Results land on the host under `~/HEAT/data/<MachFlag>_<shot>_<tag>/<timestep>/<PFCname>/`
  (VTP meshes / point clouds / CSVs, viewable in ParaView). The global log is
  `~/HEAT/data/HEATlog.txt`.

## Case directory layout

```
<case>/batchFile.dat
<case>/<MachFlag>/            # e.g. sparc/
    <equilibrium files>       # GEQDSK format
    <CAD files>               # STEP / IGES / FCStd (or STL for BYOM)
    PFCs*.csv                 # PFC definition file(s)
    <MACHINE>_input.csv       # HEAT input file
```

`MachFlag` must be one of: `sparc, arc, cmod, d3d, nstx, st40, step, west, kstar, aug, tcv, other`.

## Step-by-step: setting up a run

Prefer copying an existing case (same machine, similar physics) and modifying it over
building one from nothing. A batch file template can be generated with
`launchHEAT.py --sB <path>`.

### 1. Place the CAD

- Put STEP/IGES/FCStd files in `<case>/<MachFlag>/`. Symlinks are fine **if they are
  relative and their target is inside the case directory** (only the case dir is
  mounted into the container; absolute symlinks or targets outside it will dangle).
- Filenames with spaces survive batch-file parsing (comma-separated), but prefer
  space-free names/symlinks to avoid downstream surprises.

### 2. Verify PFC part names against the CAD

The `PFCname`, `intersectName`, and `excludeName` entries must match the **FreeCAD
import labels**, not names grepped out of the raw STEP text (FreeCAD auto-generates
labels like `COMPOUND053` for unnamed compounds). Never assume a label — list them:

```bash
docker run --rm -v <case>/<MachFlag>:/CAD --entrypoint "" plasmapotential/heat:<tag> \
  python3 -c "
import sys; sys.path.append('/usr/lib/freecad-python3/lib')
import FreeCAD, Import
Import.open('/CAD/<file>.stp')
print([o.Label for o in FreeCAD.ActiveDocument.Objects])
"
```

This takes ~1–2 min for a large (~25 MB) assembly.

### 3. Write the PFC file

CSV with header (units: `resolution` in **mm**, `timesteps` in **s**):

```
timesteps, PFCname, resolution, DivCode, intersectName, excludeName
0:10000, COMPOUND053, 0.5, LI, all, none
```

- One row per PFC to compute heat flux on; `#` comments allowed.
- `DivCode` ∈ `UI, UO, LI, LO` (upper/lower, inner/outer). It selects which `frac??`
  power fraction from the input file this PFC receives.
- `intersectName`: parts that can shadow this PFC — `all`, or a `:`-separated list of
  labels. `all` is physically safe but meshes the entire assembly (slow for big CAD).
- `excludeName`: parts to exclude from intersection checks, or `none`.

### 4. Bring the input file up to date

Input files drift as HEAT gains variables. Diff the case's input file variable names
against the canonical template `source/inputs/default_input.csv`:

```bash
diff <(grep -v '^#' source/inputs/default_input.csv | cut -d, -f1 | sed 's/ *$//' | sort) \
     <(grep -v '^#' <case>/<MachFlag>/<input>.csv | cut -d, -f1 | sed 's/ *$//' | sort)
```

Add any missing variables with the default values from the template; keep the case's
existing physics values. Sanity-check that the `fracUI/fracUO/fracLI/fracLO` values are
consistent with the `DivCode`s used in the PFC file (a PFC whose DivCode has frac=0
receives no power).

### 5. Equilibrium files

GEQDSK format, **psi in Wb/rad (already divided by 2π)**, with Bt0/Fpol/Psi/Ip signs
reflecting COCOS. The `TimeStep` column in the batch file — not the EQ filename — is
what HEAT uses for time.

### 6. Write the batch file

```
MachFlag, Tag, Shot, TimeStep, GEQDSK, CAD, PFC, Input, Output
sparc, myRun_nom, 1, 0.0, myEq.geqdsk, assy.stp, PFCs.csv, SPARC_input.csv, hfOpt:psiN:bdotn
```

- All file columns are **basenames** relative to `<case>/<MachFlag>/`.
- **Tags become directory names** — no spaces; avoid exotic characters.
- Time-varying run: repeat the tag on multiple rows, changing `TimeStep` + `GEQDSK`.
- `Output` options (`:`-separated): `hfOpt, hfGyro, hfRad, hfFil, hfRE, B, psiN, pwrDir,
  bdotn, norm, T, elmer, Btrace` — see the comment block in any existing batchFile.dat.
- The column may be named `GEQDSK` or `EQ` (both accepted).

### 7. Point docker-compose at the case

In `docker/docker-compose.yml`, set the batch-mode bind mount and ensure the `command:`
line is active:

```yaml
command: ["--m", "t", "--f", "/root/terminal/batchFile.dat"]
volumes:
  - ${HOME}/HEAT:/root/HEAT
  - /abs/path/to/<case>:/root/terminal
```

`~/HEAT` must exist on the host (output + logs). Uncomment the repo-source mount only
if the run needs unreleased code changes.

### 8. Pre-flight validation (always do this)

Parse the batch file exactly as HEAT does and check every referenced file resolves
**inside the container** (this also exercises symlinks):

```bash
docker run --rm -v <case>:/root/terminal --entrypoint "" plasmapotential/heat:<tag> \
  python3 -c "
import pandas as pd, os
os.chdir('/root/terminal')
d = pd.read_csv('batchFile.dat', sep=',', comment='#', skipinitialspace=True)
print(d.to_string())
mach = d['MachFlag'].iloc[0].strip()
eqCol = 'GEQDSK' if 'GEQDSK' in d else 'EQ'
missing = [os.path.join(mach, r[c].strip()) for _, r in d.iterrows()
           for c in [eqCol,'CAD','PFC','Input']
           if not os.path.exists(os.path.join(mach, r[c].strip()))]
print('MISSING:', missing) if missing else print('all files resolve')
"
```

### 9. Run

```bash
cd docker && docker compose up          # runs all tags in the batch file
```

or as a one-off without editing compose:

```bash
docker run --rm -v ~/HEAT:/root/HEAT -v <case>:/root/terminal \
  plasmapotential/heat:<tag> --m t --f /root/terminal/batchFile.dat
```

Watch for `Number of simulations to be scheduled from batchFile: N` early in stdout —
if N is wrong, the batch file tag/comment structure is wrong. Errors like
`Part <name> not a meshable CAD solid` or a part-not-found message mean step 2 was
skipped or wrong.

## Common pitfalls

- **PFC names are FreeCAD labels**, generated at import time — verify them (step 2).
- **Only the first row of a tag** sets CAD/PFC/Input/Output; later rows' values are
  silently ignored.
- **Stale mesh cache**: HEAT caches meshes under `~/HEAT/data/`; `overWriteMask, True`
  in the input file forces shadow-mask recompute. If CAD changed but filenames didn't,
  clear the old case output directory.
- **`intersectName: all`** on a large assembly makes meshing the dominant cost.
- **Symlinks** must stay within the mounted case directory.
- **Units**: resolution mm, timesteps s, power MW (`P` in input file), lq in mm.
- **radFile / elmerDir paths** in the input file are container paths
  (`/root/terminal/...`), not host paths.
- **`mafot_gpu, True` on a GPU MAFOT was not built for** fails at the first field-line
  trace with `CUDA error ... no kernel image is available for execution on the device`,
  then `MAFOT exited with code 1` from `MHDClass.runMAFOT`. Either set `mafot_gpu, False`
  (CPU trace; `rayTracer` may stay `mitsuba_gpu`) or rebuild the image with that GPU's
  `sm_XY` added to MAFOT's `NVCCFLAGS`.
