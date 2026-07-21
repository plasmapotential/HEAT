# docker/CLAUDE.md — maintainer notes for the HEAT docker workflow

Invariants and gotchas behind `run.sh`, `docker-compose.yml`, `.env`, and
`setup.sh`. Read this before changing any of them.

## How the pieces fit

- `../run.sh` (repo root) is the user-facing entrypoint. It validates
  configuration host-side, then runs `docker compose` **from this folder** so
  compose finds `.env` and auto-merges `docker-compose.override.yml`.
- `docker-compose.yml` parameterizes everything via `${VAR:-default}`, so the
  stack works with no `.env` at all. Keep that property: never add a variable
  without a fallback default, or bare `docker compose up` breaks for users who
  skipped `setup.sh`.
- `setup.sh` generates/updates `.env` (gitignored) from `.env.example`
  (committed). It only rewrites the machine-specific keys (`dockerUID`,
  `dockerGID`, `HEAT_DATA_DIR`), preserving user edits like a pinned tag.
- `setup.sh --dev` writes `docker-compose.override.yml` (gitignored), which
  bind-mounts the checkout over `/root/source/HEAT`. Deleting the file reverts
  to the image's baked-in source. HEAT is pure Python at runtime, so this runs
  local code changes with no image rebuild.

## ENTRYPOINT + CMD dispatch

The image defines `ENTRYPOINT launchHEAT.py` and `CMD --a 0.0.0.0 --p 8050 --m g`
(see `Dockerfile`). Consequences:

- Bare `up`/`run` → GUI. Extra args (`--m t --f ...`) replace CMD → TUI.
- A shell requires `--entrypoint ""` — `run.sh shell` does this; plain
  `docker compose run HEAT /bin/bash` would pass `/bin/bash` as launchHEAT args.
- `run.sh tui` uses `compose run --rm` (not `up`): no port publish (so it works
  while a GUI is up) and the container is removed on exit.

## Env var precedence (the subtle one)

Compose resolves a variable as: exported shell environment > `.env` > the
`:-` default in the compose file. `run.sh`'s `resolve_env()` mirrors the
first two levels only — compose-file `:-` defaults are not visible to it, so
it returns empty for a variable set nowhere else and callers supply their own
fallback (e.g. `${port:-8050}`). It deliberately uses `printenv` rather than
`$VAR`: bash sets `UID` as an unexported shell variable, so `$UID` is not
what compose (a child process) sees. Any new validation in `run.sh` must go
through `resolve_env()` or it will disagree with what compose actually does.

## The dockerUID/dockerGID contract

`source/dashGUI.py` and `source/terminalUI.py` do
`int(os.environ["dockerUID"])` to chown HEAT output. Three guards keep that
safe — keep them consistent:

1. compose defaults `dockerUID=${dockerUID:-0}` (never unset → no KeyError;
   0 = plain-docker-run root behavior).
2. `run.sh check_environment()` rejects non-numeric values (would crash
   `int()` inside the container with a murkier error).
3. `run.sh ensure_owner_env()` exports the caller's IDs when nothing set them
   (covering what the deleted `runDockerCompose` scripts did).

Escape hatch for root-owned files: `./run.sh fix_permissions`.

## The container is root — everything lives under /root

`/root/source/HEAT` (source), `/root/HEAT` (data, `heatdata` env var),
`/root/terminal` (batch files). A full non-root port (`user: ${UID}:${GID}`)
was evaluated and deferred: it requires relocating those paths in the
Dockerfile and Python. Don't add a `user:` line to compose without doing that
work; the container will fail on /root permissions.

## GPU is opt-in

The `HEAT-gpu` service (profile `gpu`) carries the NVIDIA device reservation,
which errors on hosts without the nvidia runtime (macOS, CI). `run.sh` selects
it when `HEAT_GPU=1` **and** `docker info` reports an nvidia runtime; raw
compose users run `docker compose --profile gpu up HEAT-gpu`. Keep `HEAT-gpu`
as `extends: HEAT` so the services can't drift apart.

## Where the image tag lives

`scripts/release.sh <tag>` rewrites every location and ends with a self-check
grep for stale `plasmapotential/heat:vX.Y` strings. Locations: CI
`HEAT_IMAGE_TAG` (`.github/workflows/integration-tests.yml`, which also has a
guard step asserting it matches the compose default, `.env.example`, and the
test-file defaults), the `:-` default in `docker-compose.yml`, `.env.example`,
test defaults in `tests/integrationTests/{verify,test}_nstxu_hf_rad_goldens.py`,
and the runnable examples in `tests/integrationTests/README.md` and `CLAUDE.md`.
If you add a new file with a concrete tag, add it to the `TAG_DOC_FILES` list
in `release.sh` — unless its example pairs the tag with another version string
(like `HEAT_REF=`): those examples use `<ref>`/`<tag>` placeholders instead,
because a blanket tag sed would rewrite one half of the pair and corrupt them.

## CI does not use run.sh or .env

`integration-tests.yml` runs plain `docker run -v $workspace:/root/source/HEAT`
steps (GitHub Actions can't read `.env` natively). `./run.sh test <name>`
reproduces those steps locally — keep its mount/args in sync with CI if either
changes.
