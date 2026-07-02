# HEAT integration tests

End-to-end checks that run HEAT in **terminal/batch mode** inside the published Docker image, with the repository mounted so inputs and `source/` match your working tree.

## Layout

| Path | Role |
|------|------|
| `nstxuTestCase/` | NSTX-U scenarios: `batchFile*.dat` rows plus `nstx/` inputs (GEQDSK, STEP, CSVs, `NSTXU_input*.csv`). |
| `nstxuTestCase/nstxu_hf_rad_goldens.json` | Expected stdout metrics for photon radiation (`hfRad`); used by `verify_nstxu_hf_rad_goldens.py`. |
| `D3DTestCase/` | DIII-D example case (optional in CI). |
| `ciTest.py` | Minimal smoke check that the image runs and sees the mounted repo. |
| `verify_nstxu_hf_rad_goldens.py` | Runs `batchFile_rad.dat` in Docker and compares parsed prints to the JSON goldens. |
| `test_nstxu_hf_rad_goldens.py` | Thin pytest wrapper around the verifier (needs Docker). |

Batch files list one HEAT job per active row: `MachFlag, Tag, Shot, TimeStep, EQ, CAD, PFC, Input, Output`. Paths inside CSVs are written for the container layout below.

## Container paths

CI and the verifier use:

```text
-v <repo_root>:/root/source/HEAT
```

So paths such as `/root/source/HEAT/tests/integrationTests/...` in inputs resolve to files in **your** clone. If you run HEAT in compose without mounting `source/`, you will not pick up local code changes.

## CI

GitHub Actions workflow: `.github/workflows/integration-tests.yml`.

- Image tag is set once in `env.HEAT_IMAGE_TAG` (must exist on Docker Hub as `plasmapotential/heat:<tag>`).
- Steps run `docker run … plasmapotential/heat:$HEAT_IMAGE_TAG --m t --f /root/source/HEAT/tests/integrationTests/.../batchFile….dat` except **photon radiation**, which uses the Python golden verifier (still uses Docker internally).

## Running locally

From the **repository root**, with Docker installed:

**Smoke test (mounted repo visible in container):**

```bash
docker run --rm -v "$(pwd):/root/source/HEAT" --entrypoint "" \
  plasmapotential/heat:v4.2.7 \
  python3 /root/source/HEAT/tests/integrationTests/ciTest.py
```

**Full NSTX-U batch (example: optical only):**

```bash
docker run --rm -v "$(pwd):/root/source/HEAT" \
  plasmapotential/heat:v4.2.7 \
  --m t --f /root/source/HEAT/tests/integrationTests/nstxuTestCase/batchFile_optical.dat
```

**Photon radiation with golden assertions:**

```bash
python3 tests/integrationTests/verify_nstxu_hf_rad_goldens.py \
  --workspace "$(pwd)" \
  --docker-image plasmapotential/heat:v4.2.7
```

Optional: `--log-file /tmp/hf_rad.log` saves full HEAT output. Override the image with env `HEAT_DOCKER_IMAGE` if you omit `--docker-image`.

**Same check via pytest** (host needs Python + Docker):

```bash
pytest tests/integrationTests/test_nstxu_hf_rad_goldens.py -v
```

## Updating hfRad goldens

When inputs or physics change **on purpose**, rerun the verifier (or a single batch run), read the printed “Parsed metrics”, and edit `nstxuTestCase/nstxu_hf_rad_goldens.json`. Adjust `rtol` / `atol` only if you need to allow small numerical drift across platforms.

## Adding more golden checks

Pattern to copy: a small JSON next to the case, a `verify_*.py` that `docker run`s the right `batchFile` and asserts on **stable** outputs (parsed stdout, or a small artifact under the shot path if you add one). Wire a new step into `integration-tests.yml` next to the existing jobs.
