"""Run DLUTFEM for a VGM file through an ASCII-only temporary workspace."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_SOLVER_BIN = ROOT.parent / "FEM" / "pimlfem" / "Artifact" / "Mingw-w64-gcc13.2" / "bin"
DEFAULT_WORK_ROOT = Path(r"C:\temp\dlutfem_ascii")


def is_ascii_path(path: Path) -> bool:
    try:
        str(path).encode("ascii")
    except UnicodeEncodeError:
        return False
    return True


def safe_job_name(text: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("._")
    return name or "job"


def copy_solver_runtime(source_bin: Path, target_bin: Path, refresh: bool = False) -> None:
    if refresh and target_bin.exists():
        shutil.rmtree(target_bin)
    target_bin.mkdir(parents=True, exist_ok=True)

    ignore_names = {
        "run",
        "run_tmp",
        "__pycache__",
    }
    ignore_suffixes = {
        ".db",
        ".h5",
        ".vtk",
        ".vgm",
        ".log",
        ".tmp",
    }

    for item in source_bin.iterdir():
        if item.name in ignore_names or item.suffix.lower() in ignore_suffixes:
            continue
        target = target_bin / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        elif refresh or not target.exists() or item.stat().st_mtime > target.stat().st_mtime:
            shutil.copy2(item, target)


def clear_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for item in path.iterdir():
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink()


def copy_outputs(run_dir: Path, output_dir: Path, job: str) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    for item in run_dir.iterdir():
        if item.name == "input.vgm":
            continue
        if item.is_file() and (item.name.startswith(job) or item.suffix.lower() in {".h5", ".vtk", ".db", ".log", ".txt"}):
            target = output_dir / item.name
            shutil.copy2(item, target)
            copied.append(target)
    return copied


def run_dlutfem(
    vgm_path: Path,
    solver_bin: Path,
    work_root: Path,
    output_dir: Path,
    job: str | None,
    refresh_bin: bool,
    log=print,
) -> int:
    vgm_path = vgm_path.resolve()
    solver_bin = solver_bin.resolve()
    work_root = work_root.resolve()
    output_dir = output_dir.resolve()

    if not vgm_path.exists():
        raise FileNotFoundError(f"找不到 VGM 文件: {vgm_path}")
    exe_source = solver_bin / "DLUTFEM.exe"
    if not exe_source.exists():
        raise FileNotFoundError(f"找不到 DLUTFEM.exe: {exe_source}")
    if not is_ascii_path(work_root):
        raise ValueError(f"临时运行目录必须是纯英文/ASCII 路径: {work_root}")

    job_name = safe_job_name(job or vgm_path.stem)
    runtime_bin = work_root / "bin"
    run_dir = work_root / "run" / job_name

    log(f"[1/5] 准备求解器运行库: {runtime_bin}")
    copy_solver_runtime(solver_bin, runtime_bin, refresh=refresh_bin)

    log(f"[2/5] 准备英文运行目录: {run_dir}")
    clear_directory(run_dir)
    input_vgm = run_dir / "input.vgm"
    shutil.copy2(vgm_path, input_vgm)

    exe = runtime_bin / "DLUTFEM.exe"
    cmd = [str(exe), "-i", str(input_vgm), "-w", str(run_dir), "-j", job_name]
    env = os.environ.copy()
    env["PATH"] = str(runtime_bin) + os.pathsep + env.get("PATH", "")

    log("[3/5] 开始计算:")
    log(" ".join(f'"{part}"' if " " in part else part for part in cmd))
    start = time.time()
    process = subprocess.Popen(
        cmd,
        cwd=runtime_bin,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8",
        errors="replace",
    )
    output_lines: list[str] = []
    assert process.stdout is not None
    for line in process.stdout:
        text = line.rstrip("\r\n")
        output_lines.append(text)
        log(text)
    returncode = process.wait()
    seconds = time.time() - start

    log_path = run_dir / f"{job_name}.dlutfem.log"
    log_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8", errors="replace")
    log(f"[4/5] 计算耗时: {seconds:.3f} 秒，返回码: {returncode}")

    log(f"[5/5] 回拷结果到: {output_dir}")
    copied = copy_outputs(run_dir, output_dir, job_name)
    for path in copied:
        log(f"  {path}")
    if returncode != 0:
        log(f"DLUTFEM 返回非零状态，日志已保存: {output_dir / log_path.name}")
    return returncode


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DLUTFEM for a VGM file through an ASCII-only workspace.")
    parser.add_argument("--vgm", required=True, type=Path, help="Input .vgm file.")
    parser.add_argument("--solver-bin", type=Path, default=DEFAULT_SOLVER_BIN, help="Directory containing DLUTFEM.exe.")
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK_ROOT, help="ASCII-only temporary root.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory to copy result files back to.")
    parser.add_argument("--job", default=None, help="ASCII job name. Default: sanitized VGM stem.")
    parser.add_argument("--refresh-bin", action="store_true", help="Re-copy the whole solver runtime.")
    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir is not None else args.vgm.resolve().parent / "dlutfem_result"
    raise_code = run_dlutfem(args.vgm, args.solver_bin, args.work_root, output_dir, args.job, args.refresh_bin)
    raise SystemExit(raise_code)


if __name__ == "__main__":
    main()
