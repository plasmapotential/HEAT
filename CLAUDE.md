# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is HEAT

The **Heat flux Engineering Analysis Toolkit (HEAT)** is a Python suite for predicting heat flux incident on plasma-facing components (PFCs) in tokamaks. It combines CAD geometry, MHD equilibria, and multiple heat flux models (optical, ion gyro orbit, photon radiation, filaments, runaway electrons, 3D fields) into one framework. Developed by Tom Looby at Commonwealth Fusion Systems; used to design SPARC PFCs.

## Running HEAT

HEAT runs inside Docker. The published image is `plasmapotential/heat:<tag>`; the tag lives in `docker/.env` (default in `docker/docker-compose.yml` / `docker/.env.example`) and `.github/workflows/integration-tests.yml` → `HEAT_IMAGE_TAG`, kept in sync by `scripts/release.sh`.

`./run.sh` at the repo root is the preferred entrypoint for every mode (see `docker/CLAUDE.md` for the invariants behind it):

```bash
./docker/setup.sh          # once: generate docker/.env (UID/GID, data dir, tag)
./docker/setup.sh --dev    # dev: also mount this checkout over the image's HEAT source
./run.sh gui               # GUI web app on localhost:8050 (-d for detached)
./run.sh tui path/to/batchFile.dat   # TUI/batch mode (mounts the file's dir at /root/terminal)
./run.sh shell             # interactive bash in the container
./run.sh test optical      # run an integration test case (mirrors CI)
```

Raw equivalents: `cd docker && docker compose up` (GUI); `docker compose run --rm --entrypoint "" HEAT /bin/bash` (shell). One-off overrides via exported env vars, e.g. `HEAT_IMAGE_TAG=v4.2.7 ./run.sh gui`, `HEAT_GPU=1 ./run.sh gui` (NVIDIA runtime required).

Local code changes are picked up without rebuilding the image only after `./docker/setup.sh --dev` (writes the gitignored `docker/docker-compose.override.yml` mounting this checkout at `/root/source/HEAT`); `./run.sh test` always mounts the checkout.

## Running Tests

All tests run inside the Docker container against the published image (they rely on the full compiled environment: FreeCAD, OpenFOAM, Open3D, Mitsuba).

**Integration test case via the wrapper (mounts this checkout as the HEAT source, mirrors CI):**
```bash
./run.sh test optical      # any batchFile_<name>.dat under tests/integrationTests/nstxuTestCase/
```

**Smoke test (sanity check the mount):**
```bash
docker run --rm -v "$(pwd):/root/source/HEAT" --entrypoint "" \
  plasmapotential/heat:v4.3 \
  python3 /root/source/HEAT/tests/integrationTests/ciTest.py
```

**Single integration test case, raw docker (what CI runs):**
```bash
docker run --rm -v "$(pwd):/root/source/HEAT" \
  plasmapotential/heat:v4.3 \
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
  --workspace "$(pwd)" --docker-image plasmapotential/heat:v4.3

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

This updates the image tag everywhere it lives (CI `HEAT_IMAGE_TAG`, the compose/`.env.example` defaults, test-script defaults, and doc examples), runs a self-check grep for stale tags, optionally builds the image, then prompts you to push and open a PR to `main`.

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
- **MAFOT**: external binary for field-line integration with 3D perturbed fields
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
