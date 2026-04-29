"""
Pytest wrapper for NSTX-U hfRad stdout goldens (same checks as CI).

Requires Docker and the HEAT image (HEAT_DOCKER_IMAGE, default plasmapotential/heat:v4.2.7).

Run from repo root:
  pytest tests/integrationTests/test_nstxu_hf_rad_goldens.py -v

File generated with cursor ai, reviewed by TL
"""
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VERIFY_SCRIPT = Path(__file__).resolve().parent / "verify_nstxu_hf_rad_goldens.py"


def test_nstxu_hf_rad_stdout_goldens():
    assert VERIFY_SCRIPT.is_file(), "Missing {}".format(VERIFY_SCRIPT)
    env = os.environ.copy()
    env.setdefault("HEAT_DOCKER_IMAGE", "plasmapotential/heat:v4.2.7")
    subprocess.run(
        [
            sys.executable,
            str(VERIFY_SCRIPT),
            "--workspace",
            str(REPO_ROOT),
            "--docker-image",
            env["HEAT_DOCKER_IMAGE"],
        ],
        check=True,
        env=env,
    )
