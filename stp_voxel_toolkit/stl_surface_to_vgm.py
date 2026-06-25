"""Generate VGM from an occupancy VTK and two selected surface STL files."""

from __future__ import annotations

import argparse
import math
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from vtk_selection_to_vgm import read_occupancy, write_vgm_stream


def read_ascii_or_binary_stl(path: Path) -> np.ndarray:
    data = path.read_bytes()
    if len(data) >= 84:
        tri_count = int.from_bytes(data[80:84], "little", signed=False)
        expected = 84 + tri_count * 50
        if expected == len(data):
            arr = np.frombuffer(data, dtype=np.uint8, offset=84)
            records = arr.reshape(tri_count, 50)
            vertices = np.empty((tri_count, 3, 3), dtype=np.float64)
            for i, rec in enumerate(records):
                floats = np.frombuffer(rec.tobytes()[:48], dtype="<f4")
                vertices[i] = floats[3:12].reshape(3, 3)
            return vertices

    text = data.decode("utf-8", errors="ignore")
    nums = []
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("vertex"):
            parts = stripped.split()
            if len(parts) >= 4:
                nums.append((float(parts[1]), float(parts[2]), float(parts[3])))
    if len(nums) % 3 != 0 or not nums:
        raise ValueError(f"Cannot parse STL vertices from {path}")
    return np.asarray(nums, dtype=np.float64).reshape(-1, 3, 3)


def point_triangle_distance_squared(points: np.ndarray, tri: np.ndarray) -> np.ndarray:
    # Vectorized implementation of closest point on triangle from Real-Time Collision Detection.
    a, b, c = tri
    ab = b - a
    ac = c - a
    ap = points - a
    d1 = ap @ ab
    d2 = ap @ ac
    result = np.empty(points.shape[0], dtype=np.float64)
    unresolved = np.ones(points.shape[0], dtype=bool)

    mask = (d1 <= 0.0) & (d2 <= 0.0)
    result[mask] = np.sum((points[mask] - a) ** 2, axis=1)
    unresolved &= ~mask

    bp = points - b
    d3 = bp @ ab
    d4 = bp @ ac
    mask = unresolved & (d3 >= 0.0) & (d4 <= d3)
    result[mask] = np.sum((points[mask] - b) ** 2, axis=1)
    unresolved &= ~mask

    vc = d1 * d4 - d3 * d2
    mask = unresolved & (vc <= 0.0) & (d1 >= 0.0) & (d3 <= 0.0)
    if mask.any():
        v = d1[mask] / (d1[mask] - d3[mask])
        closest = a + v[:, None] * ab
        result[mask] = np.sum((points[mask] - closest) ** 2, axis=1)
        unresolved &= ~mask

    cp = points - c
    d5 = cp @ ab
    d6 = cp @ ac
    mask = unresolved & (d6 >= 0.0) & (d5 <= d6)
    result[mask] = np.sum((points[mask] - c) ** 2, axis=1)
    unresolved &= ~mask

    vb = d5 * d2 - d1 * d6
    mask = unresolved & (vb <= 0.0) & (d2 >= 0.0) & (d6 <= 0.0)
    if mask.any():
        w = d2[mask] / (d2[mask] - d6[mask])
        closest = a + w[:, None] * ac
        result[mask] = np.sum((points[mask] - closest) ** 2, axis=1)
        unresolved &= ~mask

    va = d3 * d6 - d5 * d4
    mask = unresolved & (va <= 0.0) & ((d4 - d3) >= 0.0) & ((d5 - d6) >= 0.0)
    if mask.any():
        w = (d4[mask] - d3[mask]) / ((d4[mask] - d3[mask]) + (d5[mask] - d6[mask]))
        closest = b + w[:, None] * (c - b)
        result[mask] = np.sum((points[mask] - closest) ** 2, axis=1)
        unresolved &= ~mask

    if unresolved.any():
        denom = va[unresolved] + vb[unresolved] + vc[unresolved]
        v = vb[unresolved] / denom
        w = vc[unresolved] / denom
        closest = a + ab * v[:, None] + ac * w[:, None]
        result[unresolved] = np.sum((points[unresolved] - closest) ** 2, axis=1)

    return result


def surface_stl_to_voxels(
    stl_path: Path,
    occupancy: np.ndarray,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    tolerance: float,
    max_candidates_per_tri: int = 2_000_000,
    progress=None,
    workers: int = 1,
    chunk_size: int = 64,
) -> set[tuple[int, int, int]]:
    triangles = read_ascii_or_binary_stl(stl_path)
    nx, ny, nz = occupancy.shape
    origin_arr = np.asarray(origin, dtype=np.float64)
    spacing_arr = np.asarray(spacing, dtype=np.float64)
    tol_vec = np.array([tolerance, tolerance, tolerance], dtype=np.float64)
    tol2 = tolerance * tolerance

    def process_range(start: int, end: int) -> tuple[int, set[tuple[int, int, int]]]:
        local_selected: set[tuple[int, int, int]] = set()
        for tri_index in range(start, end):
            tri = triangles[tri_index]
            lo = tri.min(axis=0) - tol_vec
            hi = tri.max(axis=0) + tol_vec
            imin, jmin, kmin = np.floor((lo - origin_arr) / spacing_arr - 0.5).astype(int)
            imax, jmax, kmax = np.ceil((hi - origin_arr) / spacing_arr - 0.5).astype(int)
            imin, jmin, kmin = max(0, imin), max(0, jmin), max(0, kmin)
            imax, jmax, kmax = min(nx - 1, imax), min(ny - 1, jmax), min(nz - 1, kmax)
            if imin > imax or jmin > jmax or kmin > kmax:
                continue
            candidate_count = (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1)
            if candidate_count > max_candidates_per_tri:
                raise RuntimeError(
                    f"{stl_path.name} triangle {tri_index + 1} covers {candidate_count:,} voxel candidates. "
                    "Use a smaller selected-surface STL or reduce --tolerance."
                )

            ii, jj, kk = np.mgrid[imin : imax + 1, jmin : jmax + 1, kmin : kmax + 1]
            occ_mask = occupancy[ii, jj, kk].astype(bool)
            if not occ_mask.any():
                continue
            ii = ii[occ_mask].ravel()
            jj = jj[occ_mask].ravel()
            kk = kk[occ_mask].ravel()
            centers = origin_arr + (np.column_stack([ii, jj, kk]).astype(np.float64) + 0.5) * spacing_arr
            close = point_triangle_distance_squared(centers, tri) <= tol2
            for i, j, k in zip(ii[close], jj[close], kk[close]):
                local_selected.add((int(i), int(j), int(k)))
        return end - start, local_selected

    total_triangles = len(triangles)
    if total_triangles == 0:
        return set()
    workers = max(1, int(workers))
    chunk_size = max(1, int(chunk_size))
    ranges = [(start, min(start + chunk_size, total_triangles)) for start in range(0, total_triangles, chunk_size)]
    selected: set[tuple[int, int, int]] = set()
    done = 0

    if progress is not None:
        progress(0, total_triangles, 0)

    if workers == 1 or len(ranges) == 1:
        for start, end in ranges:
            count, partial = process_range(start, end)
            selected.update(partial)
            done += count
            if progress is not None:
                progress(done, total_triangles, len(selected))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(process_range, start, end) for start, end in ranges]
            for future in as_completed(futures):
                count, partial = future.result()
                selected.update(partial)
                done += count
                if progress is not None:
                    progress(done, total_triangles, len(selected))
    if progress is not None:
        progress(total_triangles, total_triangles, len(selected))
    return selected


def vector3(text: str) -> tuple[float, float, float]:
    nums = [float(x) for x in re.split(r"[, ]+", text.strip()) if x]
    if len(nums) != 3:
        raise argparse.ArgumentTypeError("expected three numbers, e.g. 0,0,-1")
    return nums[0], nums[1], nums[2]


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate VGM from occupancy VTK and selected surface STL files.")
    parser.add_argument("--occupancy", required=True, type=Path)
    parser.add_argument("--constraint-stl", required=True, type=Path)
    parser.add_argument("--load-stl", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fine", nargs=3, type=int, default=(5, 5, 5))
    parser.add_argument("--tolerance", type=float, default=None, help="Distance tolerance. Default: voxel diagonal.")
    parser.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) // 2)), help="Parallel worker threads for STL surface mapping.")
    parser.add_argument("--chunk-size", type=int, default=64, help="Triangles per parallel task.")
    parser.add_argument("--spc-dof", default="123")
    parser.add_argument("--spc-node-mode", choices=("cell-corners", "cell-min-corner"), default="cell-corners")
    parser.add_argument("--force", type=float, default=1000.0)
    parser.add_argument("--direction", type=vector3, default=(0.0, 0.0, -1.0))
    args = parser.parse_args()

    occupancy, origin, spacing = read_occupancy(args.occupancy)
    tolerance = args.tolerance
    if tolerance is None:
        tolerance = math.sqrt(spacing[0] ** 2 + spacing[1] ** 2 + spacing[2] ** 2)
    print(f"Grid: {occupancy.shape[0]}x{occupancy.shape[1]}x{occupancy.shape[2]}")
    print(f"Spacing: {spacing}, tolerance={tolerance:.9g}")

    print(f"Workers: {args.workers}, chunk_size={args.chunk_size}")
    constraint_voxels = surface_stl_to_voxels(
        args.constraint_stl,
        occupancy,
        origin,
        spacing,
        tolerance,
        workers=args.workers,
        chunk_size=args.chunk_size,
    )
    print(f"Constraint voxels: {len(constraint_voxels):,}")
    load_voxels = surface_stl_to_voxels(
        args.load_stl,
        occupancy,
        origin,
        spacing,
        tolerance,
        workers=args.workers,
        chunk_size=args.chunk_size,
    )
    print(f"Load voxels: {len(load_voxels):,}")

    def progress(stage, done, total):
        print(f"{stage}: {done:,}/{total:,}")

    write_vgm_stream(
        args.output,
        occupancy,
        origin,
        spacing,
        constraint_voxels,
        load_voxels,
        tuple(args.fine),
        args.spc_dof,
        1,
        args.force,
        args.direction,
        1,
        "Voxel VGM from selected STL surfaces",
        progress=progress,
        spc_node_mode=args.spc_node_mode,
    )
    print(f"Wrote {args.output}")
    print(f"Size: {args.output.stat().st_size / 1048576:.2f} MB")


if __name__ == "__main__":
    main()
