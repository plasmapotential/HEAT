#!/usr/bin/env python3
"""
Run NSTX-U hfRad integration (batchFile_rad.dat) in Docker and assert stdout metrics.

Usage (from repo root, with Docker available):
  python3 tests/integrationTests/verify_nstxu_hf_rad_goldens.py

CI passes HEAT image via --docker-image (see .github/workflows/integration-tests.yml).

Golden numbers live in nstxuTestCase/nstxu_hf_rad_goldens.json.

File generated with cursor ai, reviewed by TL
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path


def _repo_root_from_this_file() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_goldens_path() -> Path:
    return (
        Path(__file__).resolve().parent
        / "nstxuTestCase"
        / "nstxu_hf_rad_goldens.json"
    )


def run_heat_hf_rad(workspace: Path, docker_image: str) -> str:
    batch = "/root/source/HEAT/tests/integrationTests/nstxuTestCase/batchFile_rad.dat"
    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{workspace.resolve()}:/root/source/HEAT",
        docker_image,
        "--m",
        "t",
        "--f",
        batch,
    ]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    log = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        sys.stderr.write(log)
        raise SystemExit(
            "HEAT hfRad docker run failed with code {}.\nCommand: {}".format(
                proc.returncode, " ".join(cmd)
            )
        )
    return log


def parse_hf_rad_metrics(log: str) -> dict[str, float]:
    """Extract metrics printed by engineClass for hfRad (last match if repeated)."""
    sums = re.findall(
        r"Summation radiated power to this PFC\s*=\s*([-+eE\d.]+)", log
    )
    peaks = re.findall(r"Peak qRad to this PFC:\s*([-+eE\d.]+)", log)
    if not sums or not peaks:
        raise ValueError(
            "Could not parse hfRad summary lines from HEAT log. "
            "Expected 'Summation radiated power...' and 'Peak qRad...'."
        )

    tallies = re.findall(r"^mergedPFCs:\s+([-+eE\d.]+)", log, re.MULTILINE)
    if not tallies:
        raise ValueError(
            "Could not parse final photon tally line 'mergedPFCs:'. "
            "If the test PFC name changes, extend the parser or goldens file."
        )

    return {
        "summation_radiated_power_to_pfc": float(sums[-1]),
        "peak_qrad": float(peaks[-1]),
        "final_tally_pfc_prad_sum": float(tallies[-1]),
    }


def assert_close(
    name: str, actual: float, expected: float, rtol: float, atol: float
) -> None:
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(
            "{} out of tolerance: actual={!r} expected={!r} rtol={!r} atol={!r}".format(
                name, actual, expected, rtol, atol
            )
        )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--workspace",
        type=Path,
        default=_repo_root_from_this_file(),
        help="Repository root on host (mounted at /root/source/HEAT in container).",
    )
    p.add_argument(
        "--docker-image",
        default=os.environ.get("HEAT_DOCKER_IMAGE", "plasmapotential/heat:v4.2.7"),
        help="Docker image for HEAT (default: env HEAT_DOCKER_IMAGE or plasmapotential/heat:v4.2.7).",
    )
    p.add_argument(
        "--goldens",
        type=Path,
        default=_default_goldens_path(),
        help="JSON file with expected values and rtol.",
    )
    p.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="If set, write full HEAT stdout/stderr to this path.",
    )
    args = p.parse_args(argv)

    goldens_path = args.goldens.resolve()
    with open(goldens_path, encoding="utf-8") as f:
        spec = json.load(f)

    exp = spec["expected"]
    rtol = float(spec.get("rtol", 1e-4))
    atol = float(spec.get("atol", 0.0))

    print("Running HEAT hfRad in Docker: {}".format(args.docker_image))
    print("Workspace mount: {} -> /root/source/HEAT".format(args.workspace))
    log = run_heat_hf_rad(args.workspace, args.docker_image)
    if args.log_file:
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
        args.log_file.write_text(log, encoding="utf-8")
        print("Wrote log to {}".format(args.log_file))

    got = parse_hf_rad_metrics(log)
    print("Parsed metrics: {}".format(got))

    for key in exp:
        if key not in got:
            raise SystemExit("Golden key {!r} missing from parsed metrics".format(key))
        assert_close(key, got[key], float(exp[key]), rtol, atol)
        print("OK  {} ~ {}".format(key, exp[key]))

    print("All hfRad golden checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
