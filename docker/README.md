# HEAT Docker

## Build order and caching

The Dockerfile is ordered so that **changing the HEAT branch/tag (`HEAT_REF`) does not force a full rebuild** of OpenFOAM and swak4Foam:

1. **Stage 1 (heatbuilder)** builds, in order: base deps → M3DC1 → MAFOT → download OpenFOAM → **build OpenFOAM (Allwmake)** → **clone HEAT and build only heatFoam (wmake)** → swak4Foam.
2. The **HEAT repo is cloned only** in the step that copies `heatFoam` and runs `wmake`. So that layer is the only one that depends on `HEAT_REF`.
3. When you build with a new branch/tag (e.g. `--build-arg HEAT_REF=v4.3`), only the "clone HEAT + wmake heatFoam" layer rebuilds; the long OpenFOAM and swak4Foam layers come from cache (if you haven’t changed anything else).

Build from the **repository root**:

```bash
docker build -f docker/Dockerfile --build-arg HEAT_REF=v4.3 -t plasmapotential/heat:v4.2.6 .
```

## Using a pre-built heatbuilder (skip stage 1)

To avoid building the first stage at all (e.g. in CI or when you rarely change OpenFOAM/swak/M3DC1/MAFOT):

1. **Build and push heatbuilder occasionally** (e.g. when you bump OpenFOAM or swak4Foam):
   ```bash
   docker build -f docker/Dockerfile --target heatbuilder -t plasmapotential/heatbuilder:OF2306 .
   docker push plasmapotential/heatbuilder:OF2306
   ```
2. **In the Dockerfile**, comment out the entire first stage and start the second stage with:
   ```dockerfile
   FROM plasmapotential/heatbuilder:OF2306 AS heatbuilder
   ```
   Then the final image build only runs stage 2 (and still clones HEAT for the runtime copy). You can do this in a separate Dockerfile (e.g. `Dockerfile.quick`) or by editing the main Dockerfile when you want a fast build.

3. **Build the HEAT image** as usual; it will use the pre-built heatbuilder and only run stage 2.

This way the heavy OpenFOAM/swak4Foam build runs rarely; HEAT image builds stay fast and only depend on the current HEAT ref and stage 2.
