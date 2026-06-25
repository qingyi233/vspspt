"""Run timing comparisons for the voxelization pipeline."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
TOOL = ROOT / "stp_to_voxel_vtk.py"


def run_case(label: str, args: list[str]) -> None:
    print(f"\n## {label}")
    cmd = [sys.executable, str(TOOL), *args]
    print("Command:", " ".join(f'"{x}"' if " " in x else x for x in cmd))
    result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise SystemExit(result.returncode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(ROOT.parent / "test1.stp"))
    parser.add_argument("--resolution", default="400")
    args = parser.parse_args()

    common = ["--input", args.input, "--resolution", args.resolution, "--reuse-stl"]
    run_case("C++ voxelizer", common)
    run_case("Python voxelizer", [*common, "--python"])


if __name__ == "__main__":
    main()
