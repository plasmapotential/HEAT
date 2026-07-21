# Scripts

## release.sh — Release workflow helper

Use this script to cut a new HEAT release (Docker image tag) with less manual editing.

**What it does:**

1. Updates the Docker image tag everywhere it lives so you don’t edit multiple files by hand:
   - `.github/workflows/integration-tests.yml` (env `HEAT_IMAGE_TAG`)
   - `docker/docker-compose.yml` (the `${HEAT_IMAGE_TAG:-...}` fallback default)
   - `docker/.env.example`
   - the test defaults in `tests/integrationTests/verify_nstxu_hf_rad_goldens.py` and
     `test_nstxu_hf_rad_goldens.py`, plus the runnable examples in
     `tests/integrationTests/README.md` and `CLAUDE.md`
2. Self-checks: greps tracked files for any stale `plasmapotential/heat:vX.Y` tag and
   reports leftovers.
3. Optionally runs `docker build` with the correct `HEAT_REF` (no editing the Dockerfile).

**Typical workflow:**

1. Choose the new image tag (e.g. `v4.2.6` for minor, `v4.3.0` for major).
2. From repo root, run:
   ```bash
   ./scripts/release.sh v4.2.6 v4.3 --build
   ```
   - `v4.2.6` = Docker image tag you’ll push to Docker Hub.
   - `v4.3` = branch (or tag) to clone inside the Dockerfile (`HEAT_REF`).
   - `--build` = run the Docker build after updating files (optional).
3. If you didn’t use `--build`, run the printed `docker build` command yourself.
4. Push the image: `docker push plasmapotential/heat:<IMAGE_TAG>` (the exact command is printed by the script).
5. Commit the updated files, push your branch, open a PR to `main`. CI will use the new image tag.

**Examples:**

```bash
# Only update all refs to v4.2.6 (no build)
./scripts/release.sh v4.2.6

# Update refs and build from branch v4.3, tag image as v4.2.6
./scripts/release.sh v4.2.6 v4.3 --build

# After creating git tag v4.3.0: build from that tag and push as v4.3.0
./scripts/release.sh v4.3.0 v4.3.0 --build
```

**Prereqs:** `docker`, logged in to Docker Hub (`docker login`). Run the script **from the repository root** (the Dockerfile expects build context = repo root so that paths like `docker/buildM3DC1` resolve).

**Note:** The Dockerfile uses build-arg `HEAT_REF` (default `v4.3`), so you never need to edit it to switch branch/tag—pass it at build time. The image must be built from the **repository root** with `-f docker/Dockerfile .` (not from inside `docker/`), so that `COPY docker/buildM3DC1` etc. resolve correctly.
