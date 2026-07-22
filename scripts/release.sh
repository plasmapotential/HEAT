#!/usr/bin/env bash
# Release helper for HEAT: update image tag everywhere, then optionally build and push the Docker image.
# Usage:
#   ./scripts/release.sh <IMAGE_TAG> [HEAT_REF] [--build]
#   ./scripts/release.sh v4.2.6           # update refs only; build from default branch (v4.3)"
#   ./scripts/release.sh v4.2.6 v4.3      # update refs; HEAT_REF=v4.3 for docker build"
#   ./scripts/release.sh v4.2.6 v4.3 --build   # update refs and run docker build"
#
# Prereqs: docker, logged in to Docker Hub (docker login). Run from repo root.
# After running: commit the updated files, push your branch, open PR to main.

set -e
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [ -z "$1" ]; then
  echo "Usage: $0 <IMAGE_TAG> [HEAT_REF] [--build]"
  echo "  IMAGE_TAG  e.g. v4.2.6 or v4.3.0 (Docker image tag to push to Docker Hub)"
  echo "  HEAT_REF   optional; branch or tag to clone in Dockerfile (default: v4.3)"
  echo "  --build    run 'docker build' after updating files"
  exit 1
fi
IMAGE_TAG="$1"
HEAT_REF="v4.3"
BUILD_NOW=""
[ "$2" = "--build" ] || [ "$2" = "-b" ] && BUILD_NOW=1
[ -n "$2" ] && [ "$2" != "--build" ] && [ "$2" != "-b" ] && HEAT_REF="$2"
[ "$3" = "--build" ] || [ "$3" = "-b" ] && BUILD_NOW=1

echo "=== HEAT release helper ==="
echo "  IMAGE_TAG (Docker): $IMAGE_TAG"
echo "  HEAT_REF (git clone in Dockerfile): $HEAT_REF"
echo ""

# 1) Update integration-tests workflow (single env var at top)
CI_FILE=".github/workflows/integration-tests.yml"
if [ -f "$CI_FILE" ]; then
  if grep -q "HEAT_IMAGE_TAG:" "$CI_FILE"; then
    if sed --version 2>/dev/null | grep -q GNU; then
      sed -i "s/HEAT_IMAGE_TAG: .*/HEAT_IMAGE_TAG: $IMAGE_TAG/" "$CI_FILE"
    else
      sed -i.bak "s/HEAT_IMAGE_TAG: .*/HEAT_IMAGE_TAG: $IMAGE_TAG/" "$CI_FILE" && rm -f "${CI_FILE}.bak"
    fi
    echo "Updated $CI_FILE -> HEAT_IMAGE_TAG: $IMAGE_TAG"
  else
    echo "Warning: $CI_FILE has no HEAT_IMAGE_TAG env; add it manually."
  fi
else
  echo "Warning: $CI_FILE not found."
fi

# 2) Update docker-compose image tags
for COMPOSE in docker/docker-compose.yml docker/docker-compose_local.yml docker/docker-compose_tokens.yml; do
  if [ -f "$COMPOSE" ]; then
    if grep -q 'plasmapotential/heat:' "$COMPOSE"; then
      if sed --version 2>/dev/null | grep -q GNU; then
        sed -i "s|plasmapotential/heat:[^\"']*|plasmapotential/heat:$IMAGE_TAG|g" "$COMPOSE"
      else
        sed -i.bak "s|plasmapotential/heat:[^\"']*|plasmapotential/heat:$IMAGE_TAG|g" "$COMPOSE" && rm -f "${COMPOSE}.bak"
      fi
      echo "Updated $COMPOSE -> image tag $IMAGE_TAG"
    fi
  fi
done

echo ""
echo "=== Next steps ==="
echo "1. Build the Docker image (must run from repo root; context is ., Dockerfile uses docker/ for COPYs):"
echo ""
echo "   docker build --build-arg HEAT_REF=$HEAT_REF -t plasmapotential/heat:$IMAGE_TAG -f docker/Dockerfile ."
echo ""
echo "2. Tag as latest (optional):"
echo "   docker tag plasmapotential/heat:$IMAGE_TAG plasmapotential/heat:latest"
echo ""
echo "3. Push to Docker Hub:"
echo "   docker push plasmapotential/heat:$IMAGE_TAG"
echo "   # docker push plasmapotential/heat:latest  # if you tagged latest"
echo ""
echo "4. Commit the updated CI and compose files, push your branch, and open a PR to main."
echo "   CI will run using image $IMAGE_TAG (ensure it exists on Docker Hub before merging)."
echo ""

if [ -n "$BUILD_NOW" ]; then
  echo "Running docker build (--build passed)..."
  docker build --build-arg HEAT_REF="$HEAT_REF" -t "plasmapotential/heat:$IMAGE_TAG" -f docker/Dockerfile .
  echo ""
  echo "Build done. Push with: docker push plasmapotential/heat:$IMAGE_TAG"
fi
