#!/bin/bash
# Helper for running the HEAT docker compose services.
#
# Usage:
#   Source this file and call run_heat, or run it directly:
#     ./run.sh <command> [arguments]
#
# Examples:
#   ./run.sh gui                   # start the GUI on http://localhost:8050
#   ./run.sh gui -d                # same, detached
#   ./run.sh tui path/to/batchFile.dat   # terminal/batch mode
#   ./run.sh shell                 # interactive shell in the container
#   ./run.sh test optical          # run an integration test case locally
#   ./run.sh fix_permissions       # chown the HEAT data dir back to your user
#   ./run.sh cleanup               # tidy up orphaned/stopped containers
#   ./run.sh rm_docker             # remove all images, containers and volumes
#
# Configuration lives in docker/.env (run docker/setup.sh to generate it);
# an exported shell variable overrides .env for one-off runs, e.g.:
#   HEAT_IMAGE_TAG=v4.2.7 ./run.sh gui
#   HEAT_GPU=1 ./run.sh tui batchFile.dat   # enable the NVIDIA GPU reservation

# Define color codes for pretty output
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
NOCOLOR='\033[0m' # Reset to no color

print_message() {
    local color="$1"
    local message="$2"
    printf "${color}%s${NOCOLOR}\n" "$message"
}

warn()   { print_message "${LIGHTRED}" "$1"; }
notice() { print_message "${LIGHTGREEN}" "$1"; }

HEAT_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HEAT_DOCKER_DIR="${HEAT_REPO_ROOT}/docker"

# Commands that launch the HEAT container
VALID_COMMANDS="gui tui shell test fix_permissions"
# Commands that require at least one argument
ARG_REQUIRED_COMMANDS="tui"

# Return 0 if $1 appears as a whitespace-separated word in $2
contains_word() {
  case " $2 " in
    *" $1 "*) return 0 ;;
    *) return 1 ;;
  esac
}

# All compose invocations run from the docker/ folder so that compose finds
# .env and auto-merges docker-compose.override.yml (dev source mount).
heat_compose() {
  (cd "$HEAT_DOCKER_DIR" && docker compose "$@")
}

# Resolve a variable the way docker compose does: a value exported in the
# environment wins, otherwise fall back to the value assigned in docker/.env.
# Uses printenv (not $VAR) deliberately -- bash always sets UID/EUID as shell
# variables but does not export them, so they are NOT what compose (a child
# process) actually sees; printenv reports only the exported environment.
resolve_env() {
  local key="$1" value
  value=$(printenv "$key" 2>/dev/null)
  if [ -n "$value" ]; then
    printf '%s' "$value"
    return 0
  fi
  sed -n "s/^${key}=//p" "${HEAT_DOCKER_DIR}/.env" 2>/dev/null | tail -n1
}

run_heat_help() {
  notice "Usage: run_heat <command> [arguments]"
  echo
  echo "Commands:"
  echo "  gui [-d]         Start the HEAT GUI (Dash web app on http://localhost:8050)."
  echo "                   -d runs it detached."
  echo "  tui <batchFile>  Run HEAT in terminal/batch mode. The batch file's folder is"
  echo "                   mounted at /root/terminal inside the container. Set"
  echo "                   HEAT_RUNS_DIR in docker/.env to also mount a shared runs"
  echo "                   folder at /root/HEAT_runs (for absolute input paths)."
  echo "  shell            Interactive shell in the container (bash, not launchHEAT)."
  echo "  test [name]      Run an integration test case with this checkout mounted as"
  echo "                   the HEAT source (mirrors CI). Default: optical."
  echo "  fix_permissions  chown the HEAT data dir back to your user (root escape hatch)."
  echo
  echo "Maintenance:"
  echo "  cleanup          Stop and remove orphaned/stopped containers from this project."
  echo "  rm_docker        Remove ALL containers, volumes and images for this project."
  echo
  echo "Configuration: docker/.env (generate with docker/setup.sh; setup.sh --dev also"
  echo "mounts this checkout over the image's HEAT source). An exported variable"
  echo "overrides .env, e.g. HEAT_IMAGE_TAG=v4.2.7 ./run.sh gui, HEAT_GPU=1 ./run.sh gui."
  echo
  echo "Deprecated aliases: runDockerCompose (= gui), runDockerComposeDetached"
  echo "(= detached shell container), kept for compatibility with the old scripts."
}

# Make sure the environment is set up correctly before running anything.
# Returns non-zero (and explains) if something is missing.
check_environment() {
  if [ ! -f "${HEAT_DOCKER_DIR}/docker-compose.yml" ]; then
    warn "Error: no docker-compose.yml in ${HEAT_DOCKER_DIR}."
    warn "Run this from a HEAT checkout (run.sh lives at the repo root)."
    return 1
  fi
  if [ ! -f "${HEAT_DOCKER_DIR}/.env" ]; then
    notice "Note: no docker/.env found; using compose defaults."
    notice "Run 'docker/setup.sh' to generate one (records your UID/GID so HEAT"
    notice "output is owned by you, and pins the data dir and image tag)."
  fi
  # dashGUI.py/terminalUI.py do int(os.environ["dockerUID"]) inside the
  # container, so a non-numeric value crashes HEAT at startup with a murkier
  # error. Empty is fine (compose substitutes 0 = current root behavior).
  local uid gid bad=""
  uid=$(resolve_env dockerUID)
  gid=$(resolve_env dockerGID)
  case "$uid" in *[!0-9]*) bad="dockerUID='${uid}'" ;; esac
  case "$gid" in *[!0-9]*) bad="${bad:+$bad, }dockerGID='${gid}'" ;; esac
  if [ -n "$bad" ]; then
    warn "Error: ${bad} is not a non-negative integer (set in docker/.env or the environment)."
    warn "HEAT reads dockerUID/dockerGID inside the container to chown its output."
    warn "Run 'docker/setup.sh' to regenerate docker/.env."
    return 1
  fi
  return 0
}

# Data dir: compose bind-mounts it at /root/HEAT; create it if missing so
# docker does not create it root-owned. Called only for commands that launch
# a container (not for maintenance commands or invalid input).
ensure_data_dir() {
  local data_dir
  data_dir=$(resolve_env HEAT_DATA_DIR)
  data_dir="${data_dir:-$HOME/HEAT}"
  if [ ! -d "$data_dir" ]; then
    notice "Creating HEAT data directory: ${data_dir}"
    mkdir -p "$data_dir" || { warn "Error: could not create ${data_dir}."; return 1; }
  fi
  return 0
}

# Same treatment for the optional HEAT_RUNS_DIR mount (/root/HEAT_runs).
# Empty means the compose fallback (an empty named volume) — nothing to create.
ensure_runs_dir() {
  local runs_dir
  runs_dir=$(resolve_env HEAT_RUNS_DIR)
  if [ -n "$runs_dir" ] && [ ! -d "$runs_dir" ]; then
    notice "Creating HEAT runs directory: ${runs_dir}"
    mkdir -p "$runs_dir" || { warn "Error: could not create ${runs_dir}."; return 1; }
  fi
  return 0
}

# The old runDockerCompose scripts exported dockerUID/dockerGID so HEAT chowns
# its output to the host user instead of root. Preserve that: if neither the
# environment nor .env provides them, export the caller's IDs (docker group
# GID when available, matching the old script).
ensure_owner_env() {
  if [ -z "$(resolve_env dockerUID)" ]; then
    dockerUID="$(id -u)"
    export dockerUID
  fi
  if [ -z "$(resolve_env dockerGID)" ]; then
    if command -v getent >/dev/null 2>&1 && getent group docker >/dev/null 2>&1; then
      dockerGID="$(getent group docker | cut -d: -f3)"
    else
      dockerGID="$(id -g)"
    fi
    export dockerGID
  fi
}

# Select the compose service: HEAT-gpu (NVIDIA reservation) when HEAT_GPU=1,
# plain HEAT otherwise. A GPU that was requested but is unavailable is an
# error (a silent CPU fallback would waste a long run); a GPU that is
# available but unrequested gets a warning so it isn't left idle unnoticed.
# No --profile flag is needed: compose auto-enables a service's profiles when
# it is named explicitly, and every call site below names $HEAT_SERVICE.
select_service() {
  HEAT_SERVICE="HEAT"
  local gpu_requested runtimes
  gpu_requested="$(resolve_env HEAT_GPU)"
  runtimes=$(docker info --format '{{json .Runtimes}}' 2>/dev/null)
  if [ -z "$runtimes" ]; then
    # docker info unavailable (daemon down?); skip the runtime check and let
    # compose report the real problem.
    [ "$gpu_requested" = "1" ] && HEAT_SERVICE="HEAT-gpu"
    return 0
  fi
  if [ "$gpu_requested" = "1" ]; then
    if printf '%s' "$runtimes" | grep -q nvidia; then
      HEAT_SERVICE="HEAT-gpu"
    else
      warn "Error: HEAT_GPU=1 but 'docker info' reports no nvidia runtime."
      warn "Install the NVIDIA container toolkit, or unset HEAT_GPU (in docker/.env"
      warn "or the environment) to run without a GPU."
      return 1
    fi
  elif printf '%s' "$runtimes" | grep -q nvidia; then
    warn "Note: an NVIDIA runtime is available but HEAT_GPU is not enabled; running WITHOUT the GPU."
    warn "Set HEAT_GPU=1 (in docker/.env or the environment) to use it."
  fi
  return 0
}

# Stop and remove orphaned/stopped containers left over from this project.
docker_cleanup() {
  notice "Stopping this project's containers and removing orphans..."
  heat_compose down --remove-orphans
  notice "Removing any remaining stopped 'run' containers for this project..."
  heat_compose rm --force --stop 2>/dev/null
  notice "Cleanup complete."
}

# Remove everything associated with this project's docker images:
# containers, named/anonymous volumes and the images themselves.
docker_rm() {
  # Image repository for this project (all tags of it are removed).
  local repo
  repo=$(resolve_env HEAT_IMAGE)
  repo="${repo:-plasmapotential/heat}"

  warn "This will remove ALL containers, volumes and images (every tag) for this project:"
  echo "  - ${repo}"
  printf "Are you sure? [y/N] "
  read -r reply
  case "$reply" in
    [Yy]*) ;;
    *) notice "Aborted."; return 0 ;;
  esac

  notice "Stopping containers and removing volumes/orphans..."
  heat_compose down --volumes --remove-orphans --rmi all
  notice "Removing any remaining stopped 'run' containers for this project..."
  heat_compose rm --force --stop --volumes 2>/dev/null

  # Remove every tag of the image, in case they weren't created via compose.
  notice "Removing images (all tags)..."
  local images
  images=$(docker images --filter "reference=${repo}" --quiet | sort -u)
  if [ -n "$images" ]; then
    # shellcheck disable=SC2086
    docker image rm --force $images 2>/dev/null
  fi

  notice "Purge complete."
  echo
  echo "Note: this removes only this project's images/containers/volumes."
  echo "To reclaim build cache and everything else Docker has accumulated"
  echo "(affects ALL Docker projects on this machine), run:"
  echo "  docker system prune -a --volumes"
}

run_heat() {
  local command="$1"

  if [ -z "$command" ] || [ "$command" = "help" ] || [ "$command" = "-h" ] || [ "$command" = "--help" ]; then
    run_heat_help
    return 0
  fi

  # Deprecated aliases covering the deleted docker/runDockerCompose* scripts
  if [ "$command" = "runDockerCompose" ]; then
    warn "'runDockerCompose' is deprecated; use './run.sh gui' instead."
    command="gui"
  elif [ "$command" = "runDockerComposeDetached" ]; then
    warn "'runDockerComposeDetached' is deprecated; use './run.sh shell' (interactive)"
    warn "or './run.sh gui -d' (detached GUI) instead."
    check_environment || return 1
    ensure_data_dir || return 1
    ensure_runs_dir || return 1
    ensure_owner_env
    notice "Starting a detached shell container (old runDockerComposeDetached behavior)..."
    heat_compose run -d --entrypoint /bin/bash HEAT -c "sleep infinity"
    return $?
  fi

  # Maintenance shortcuts
  if [ "$command" = "cleanup" ]; then
    check_environment || return 1
    docker_cleanup
    return $?
  fi
  if [ "$command" = "rm_docker" ]; then
    check_environment || return 1
    docker_rm
    return $?
  fi

  # Make sure the environment is set up correctly
  check_environment || return 1

  if ! contains_word "$command" "$VALID_COMMANDS"; then
    warn "Error: '$command' is not a valid command."
    echo
    run_heat_help
    return 1
  fi

  # Drop the command name so the rest are its arguments
  shift

  if contains_word "$command" "$ARG_REQUIRED_COMMANDS" && [ "$#" -eq 0 ]; then
    warn "Error: the '$command' command requires an argument."
    warn "For example: run_heat tui path/to/batchFile.dat"
    return 1
  fi

  ensure_data_dir || return 1
  ensure_runs_dir || return 1
  ensure_owner_env
  select_service || return 1

  case "$command" in
    gui)
      local detach=""
      if [ "$1" = "-d" ] || [ "$1" = "--detach" ]; then
        detach="-d"
      elif [ -n "$1" ]; then
        warn "Error: unknown gui option '$1' (only -d/--detach is supported)."
        return 1
      fi
      local port
      port=$(resolve_env HEAT_PORT)
      if [ -n "$detach" ]; then
        notice "Starting the HEAT GUI detached. Open http://localhost:${port:-8050}"
        notice "Stop it with './run.sh cleanup' (or 'docker compose down' from docker/)."
      else
        notice "Starting the HEAT GUI (Ctrl-C to stop). Open http://localhost:${port:-8050}"
      fi
      # shellcheck disable=SC2086
      heat_compose up $detach "$HEAT_SERVICE"
      ;;
    tui)
      local batch="$1"
      if [ ! -f "$batch" ]; then
        warn "Error: batch file '$batch' not found."
        return 1
      fi
      local batch_dir base
      batch_dir="$(cd "$(dirname "$batch")" && pwd)"
      base="$(basename "$batch")"
      # The batch file's folder is mounted at /root/terminal, so input files
      # referenced with /root/terminal/... paths sit alongside the batch file.
      export BATCH_DIR="$batch_dir"
      notice "Running HEAT in terminal mode on ${batch_dir}/${base}"
      heat_compose run --rm "$HEAT_SERVICE" --m t --f "/root/terminal/${base}"
      ;;
    shell)
      # Override the launchHEAT.py entrypoint so we get bash
      heat_compose run --rm --entrypoint "" "$HEAT_SERVICE" /bin/bash
      ;;
    test)
      local name="${1:-optical}"
      local case_dir="${HEAT_REPO_ROOT}/tests/integrationTests/nstxuTestCase"
      local batch="${case_dir}/batchFile_${name}.dat"
      if [ ! -f "$batch" ]; then
        warn "Error: no such test case '${name}'. Available:"
        for f in "${case_dir}"/batchFile_*.dat; do
          f="$(basename "$f")"; f="${f#batchFile_}"; echo "  ${f%.dat}"
        done
        return 1
      fi
      # Mount this checkout as the HEAT source, mirroring CI (integration-tests.yml)
      notice "Running integration test '${name}' with ${HEAT_REPO_ROOT} mounted as /root/source/HEAT"
      heat_compose run --rm -v "${HEAT_REPO_ROOT}:/root/source/HEAT" "$HEAT_SERVICE" \
        --m t --f "/root/source/HEAT/tests/integrationTests/nstxuTestCase/batchFile_${name}.dat"
      ;;
    fix_permissions)
      local uid gid
      uid=$(resolve_env dockerUID)
      gid=$(resolve_env dockerGID)
      notice "chowning /root/HEAT (the mounted HEAT data dir) to ${uid}:${gid}..."
      heat_compose run --rm --entrypoint "" HEAT chown -R "${uid}:${gid}" /root/HEAT
      ;;
  esac
}

# If executed directly (not sourced), run the function with the given arguments.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  run_heat "$@"
fi
