"""Minimal STEP/STP to voxel VTK converter.

Outputs:
1. Structured occupancy VTK: efficient for high-resolution voxel data.
2. Optional visible surface VTK: directly visible in ParaView.

Example:
    python stp_to_voxel_vtk.py --input ..\\test1.stp --target-cells 100000000 --surface
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np


ROOT = Path(__file__).resolve().parent
RUNTIME_DIR = ROOT / "pyvoxel_runtime"
if str(RUNTIME_DIR) not in sys.path:
    sys.path.insert(0, str(RUNTIME_DIR))

import gmsh
import pvoxel
import voxel


@dataclass
class Stage:
    name: str
    seconds: float
    detail: str


class StepMeshingError(RuntimeError):
    pass


def parse_positive_int(text: str) -> int:
    cleaned = str(text).strip().replace(",", "")
    value = float(cleaned)
    if not value.is_integer():
        raise argparse.ArgumentTypeError(f"expected an integer value, got {text!r}")
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError(f"expected a positive integer, got {text!r}")
    return result


def write_ascii_stl(stl_path: Path, name: str, vertices: list[tuple[float, float, float]]) -> None:
    def normal(a, b, c):
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        length = (nx * nx + ny * ny + nz * nz) ** 0.5
        if length == 0:
            return 0.0, 0.0, 0.0
        return nx / length, ny / length, nz / length

    with stl_path.open("w", encoding="ascii", newline="\n") as f:
        f.write(f"solid {name}\n")
        for i in range(0, len(vertices), 3):
            tri = vertices[i : i + 3]
            nx, ny, nz = normal(tri[0], tri[1], tri[2])
            f.write(f"  facet normal {nx:.9g} {ny:.9g} {nz:.9g}\n")
            f.write("    outer loop\n")
            for x, y, z in tri:
                f.write(f"      vertex {x:.9g} {y:.9g} {z:.9g}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {name}\n")


def collect_gmsh_triangles() -> list[tuple[float, float, float]]:
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    node_lookup = {
        int(tag): (float(coords[3 * i]), float(coords[3 * i + 1]), float(coords[3 * i + 2]))
        for i, tag in enumerate(node_tags)
    }

    vertices: list[tuple[float, float, float]] = []
    element_types, _element_tags, element_node_tags = gmsh.model.mesh.getElements(dim=2)
    for element_type, node_tags_for_type in zip(element_types, element_node_tags):
        _name, _dim, _order, node_count, _local_coords, _primary = gmsh.model.mesh.getElementProperties(element_type)
        nodes = [int(x) for x in node_tags_for_type]
        if node_count == 3:
            for i in range(0, len(nodes), 3):
                vertices.extend((node_lookup[nodes[i]], node_lookup[nodes[i + 1]], node_lookup[nodes[i + 2]]))
        elif node_count == 4:
            for i in range(0, len(nodes), 4):
                a, b, c, d = [node_lookup[nodes[i + j]] for j in range(4)]
                vertices.extend((a, b, c, a, c, d))
        else:
            raise RuntimeError(f"Unsupported surface element type {element_type} with {node_count} nodes")
    return vertices


def step_to_stl(step_path: Path, stl_path: Path, mesh_size: float | None) -> int:
    with TemporaryDirectory(prefix="stp_voxel_") as temp_name:
        temp_dir = Path(temp_name)
        temp_step = temp_dir / step_path.name
        temp_stl = temp_dir / stl_path.name
        shutil.copy2(step_path, temp_step)

        gmsh.initialize(interruptible=False)
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add(temp_step.stem)
            imported = gmsh.model.occ.importShapes(str(temp_step))
            if not imported:
                raise RuntimeError(f"Gmsh did not import any shape from {step_path}")
            gmsh.model.occ.synchronize()

            if mesh_size is not None:
                gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size)
                gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

            gmsh.option.setNumber("Mesh.Binary", 0)
            gmsh.option.setNumber("Mesh.SaveAll", 1)
            try:
                gmsh.model.mesh.generate(2)
            except Exception as exc:
                raise StepMeshingError(
                    "STEP -> STL failed: Gmsh could not mesh the CAD surfaces.\n"
                    f"Original error: {exc}\n\n"
                    "This is usually caused by problematic CAD topology such as periodic surfaces, bad faces, "
                    "degenerated edges, or unsewn faces. Export STL from CAD software and use the STL as input, "
                    "or repair/simplify the STEP model first."
                ) from exc
            triangles = collect_gmsh_triangles()
            write_ascii_stl(temp_stl, temp_step.stem, triangles)
        finally:
            gmsh.finalize()

        if not temp_stl.exists() or temp_stl.stat().st_size == 0:
            raise RuntimeError("STL writer produced an empty file")
        shutil.copy2(temp_stl, stl_path)
        return len(triangles) // 3


def grid_for_mesh(mesh_vertices: np.ndarray, target_cells: int, longest_resolution: int | None) -> tuple[int, int, int]:
    spans = np.array([mesh_vertices[:, axis, :].max() - mesh_vertices[:, axis, :].min() for axis in range(3)], dtype=float)
    max_span = float(spans.max())
    if longest_resolution is not None:
        dims = np.maximum(2, np.rint(longest_resolution * spans / max_span).astype(int))
    else:
        ratio = float(np.prod(spans / max_span))
        base = (target_cells / ratio) ** (1.0 / 3.0)
        dims = np.maximum(2, np.rint(base * spans / max_span).astype(int))
    return int(dims[0]), int(dims[1]), int(dims[2])


def write_occupancy_vtk(path: Path, occupancy: np.ndarray, origin, spacing) -> None:
    nx, ny, nz = occupancy.shape
    data = np.ascontiguousarray(occupancy.astype(np.uint8)).ravel(order="F")
    header = (
        "# vtk DataFile Version 3.0\n"
        "High-resolution voxel occupancy\n"
        "BINARY\n"
        "DATASET STRUCTURED_POINTS\n"
        f"DIMENSIONS {nx + 1} {ny + 1} {nz + 1}\n"
        f"ORIGIN {origin[0]:.17g} {origin[1]:.17g} {origin[2]:.17g}\n"
        f"SPACING {spacing[0]:.17g} {spacing[1]:.17g} {spacing[2]:.17g}\n"
        f"CELL_DATA {nx * ny * nz}\n"
        "SCALARS occupancy unsigned_char 1\n"
        "LOOKUP_TABLE default\n"
    )
    with path.open("wb") as f:
        f.write(header.encode("ascii"))
        f.write(data.tobytes())


def boundary_masks(grid: np.ndarray) -> list[tuple[str, np.ndarray]]:
    occ = grid.astype(bool, copy=False)
    masks = []
    m = occ.copy(); m[1:, :, :] &= ~occ[:-1, :, :]; masks.append(("xmin", m))
    m = occ.copy(); m[:-1, :, :] &= ~occ[1:, :, :]; masks.append(("xmax", m))
    m = occ.copy(); m[:, 1:, :] &= ~occ[:, :-1, :]; masks.append(("ymin", m))
    m = occ.copy(); m[:, :-1, :] &= ~occ[:, 1:, :]; masks.append(("ymax", m))
    m = occ.copy(); m[:, :, 1:] &= ~occ[:, :, :-1]; masks.append(("zmin", m))
    m = occ.copy(); m[:, :, :-1] &= ~occ[:, :, 1:]; masks.append(("zmax", m))
    return masks


def face_points(direction: str, ijk, origin, spacing) -> np.ndarray:
    i, j, k = ijk
    x0, y0, z0 = origin
    dx, dy, dz = spacing
    xi = x0 + i.astype(np.float64) * dx
    xip = xi + dx
    yj = y0 + j.astype(np.float64) * dy
    yjp = yj + dy
    zk = z0 + k.astype(np.float64) * dz
    zkp = zk + dz
    pts = np.empty((i.size, 4, 3), dtype=">f4")

    if direction == "xmin":
        pts[:, 0, :] = np.column_stack([xi, yj, zk]); pts[:, 1, :] = np.column_stack([xi, yj, zkp])
        pts[:, 2, :] = np.column_stack([xi, yjp, zkp]); pts[:, 3, :] = np.column_stack([xi, yjp, zk])
    elif direction == "xmax":
        pts[:, 0, :] = np.column_stack([xip, yj, zk]); pts[:, 1, :] = np.column_stack([xip, yjp, zk])
        pts[:, 2, :] = np.column_stack([xip, yjp, zkp]); pts[:, 3, :] = np.column_stack([xip, yj, zkp])
    elif direction == "ymin":
        pts[:, 0, :] = np.column_stack([xi, yj, zk]); pts[:, 1, :] = np.column_stack([xip, yj, zk])
        pts[:, 2, :] = np.column_stack([xip, yj, zkp]); pts[:, 3, :] = np.column_stack([xi, yj, zkp])
    elif direction == "ymax":
        pts[:, 0, :] = np.column_stack([xi, yjp, zk]); pts[:, 1, :] = np.column_stack([xi, yjp, zkp])
        pts[:, 2, :] = np.column_stack([xip, yjp, zkp]); pts[:, 3, :] = np.column_stack([xip, yjp, zk])
    elif direction == "zmin":
        pts[:, 0, :] = np.column_stack([xi, yj, zk]); pts[:, 1, :] = np.column_stack([xi, yjp, zk])
        pts[:, 2, :] = np.column_stack([xip, yjp, zk]); pts[:, 3, :] = np.column_stack([xip, yj, zk])
    elif direction == "zmax":
        pts[:, 0, :] = np.column_stack([xi, yj, zkp]); pts[:, 1, :] = np.column_stack([xip, yj, zkp])
        pts[:, 2, :] = np.column_stack([xip, yjp, zkp]); pts[:, 3, :] = np.column_stack([xi, yjp, zkp])
    else:
        raise ValueError(direction)
    return pts.reshape(-1, 3)


def write_surface_vtk(path: Path, grid: np.ndarray, origin, spacing) -> tuple[int, int]:
    masks = boundary_masks(grid)
    face_count = sum(int(mask.sum()) for _, mask in masks)
    point_count = face_count * 4
    nx, ny, nz = grid.shape

    def write_int_scalar(f, name: str, values: np.ndarray) -> None:
        f.write(f"SCALARS {name} int 1\nLOOKUP_TABLE default\n".encode("ascii"))
        f.write(np.asarray(values, dtype=">i4").tobytes())
        f.write(b"\n")

    with path.open("wb") as f:
        f.write((
            "# vtk DataFile Version 3.0\n"
            "Occupied voxel exterior surface\n"
            "BINARY\n"
            "DATASET UNSTRUCTURED_GRID\n"
            f"POINTS {point_count} float\n"
        ).encode("ascii"))
        for direction, mask in masks:
            if mask.any():
                f.write(face_points(direction, np.nonzero(mask), origin, spacing).tobytes())
        f.write(b"\n")
        f.write(f"CELLS {face_count} {face_count * 5}\n".encode("ascii"))
        offset = 0
        remaining = face_count
        while remaining:
            n = min(1_000_000, remaining)
            ids = np.arange(offset, offset + n * 4, dtype=np.int32).reshape(n, 4)
            polys = np.empty((n, 5), dtype=">i4")
            polys[:, 0] = 4
            polys[:, 1:] = ids.astype(">i4", copy=False)
            f.write(polys.tobytes())
            offset += n * 4
            remaining -= n
        f.write(b"\n")
        f.write(f"CELL_TYPES {face_count}\n".encode("ascii"))
        remaining = face_count
        while remaining:
            n = min(1_000_000, remaining)
            f.write(np.full(n, 9, dtype=">i4").tobytes())
            remaining -= n
        f.write(b"\n")

        voxel_i_parts = []
        voxel_j_parts = []
        voxel_k_parts = []
        voxel_eid_parts = []
        face_id_parts = []
        for face_id, (_direction, mask) in enumerate(masks):
            if not mask.any():
                continue
            i, j, k = np.nonzero(mask)
            voxel_i_parts.append(i.astype(np.int32, copy=False))
            voxel_j_parts.append(j.astype(np.int32, copy=False))
            voxel_k_parts.append(k.astype(np.int32, copy=False))
            voxel_eid_parts.append((1 + i + nx * j + nx * ny * k).astype(np.int32, copy=False))
            face_id_parts.append(np.full(i.size, face_id, dtype=np.int32))

        f.write(f"CELL_DATA {face_count}\n".encode("ascii"))
        write_int_scalar(f, "voxel_i", np.concatenate(voxel_i_parts) if voxel_i_parts else np.empty(0, dtype=np.int32))
        write_int_scalar(f, "voxel_j", np.concatenate(voxel_j_parts) if voxel_j_parts else np.empty(0, dtype=np.int32))
        write_int_scalar(f, "voxel_k", np.concatenate(voxel_k_parts) if voxel_k_parts else np.empty(0, dtype=np.int32))
        write_int_scalar(f, "voxel_eid", np.concatenate(voxel_eid_parts) if voxel_eid_parts else np.empty(0, dtype=np.int32))
        write_int_scalar(f, "face_id", np.concatenate(face_id_parts) if face_id_parts else np.empty(0, dtype=np.int32))
    return face_count, point_count


def run(args: argparse.Namespace) -> list[Stage]:
    stages: list[Stage] = []
    total_start = time.perf_counter()
    input_path = args.input.resolve()
    out_dir = args.output.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    stl_path = out_dir / f"{input_path.stem}.stl"
    if input_path.suffix.lower() in {".stp", ".step"}:
        start = time.perf_counter()
        if args.reuse_stl and stl_path.exists() and stl_path.stat().st_size > 0:
            stages.append(Stage("STEP -> STL", 0.0, f"reused {stl_path.name}"))
        else:
            tri_count = step_to_stl(input_path, stl_path, args.mesh_size)
            stages.append(Stage("STEP -> STL", time.perf_counter() - start, f"{tri_count:,} triangles, {stl_path.stat().st_size / 1048576:.2f} MB"))
    elif input_path.suffix.lower() == ".stl":
        stl_path = input_path
        stages.append(Stage("STEP -> STL", 0.0, "input is STL"))
    else:
        raise ValueError("Input must be .stp, .step, or .stl")

    start = time.perf_counter()
    mesh_vertices = pvoxel.Read_stl(str(stl_path))[0]
    mins = tuple(float(mesh_vertices[:, axis, :].min()) for axis in range(3))
    maxs = tuple(float(mesh_vertices[:, axis, :].max()) for axis in range(3))
    spans = tuple(maxs[i] - mins[i] for i in range(3))
    nx, ny, nz = grid_for_mesh(mesh_vertices, args.target_cells, args.resolution)
    total_cells = nx * ny * nz
    stages.append(Stage("Read STL / Grid", time.perf_counter() - start, f"{mesh_vertices.shape[0]:,} triangles, grid {nx}x{ny}x{nz}={total_cells:,}"))

    start = time.perf_counter()
    if args.python:
        vox_obj = pvoxel.Voxel(nx, ny, nz)
        grid = vox_obj.VOXELISE(mesh_vertices, ray=args.ray, parallel=False)
        backend = "Python"
    else:
        vox_obj = voxel.Voxel(nx, ny, nz)
        grid = vox_obj.VOXELISE(mesh_vertices, ray=args.ray, parallel=False)
        backend = "C++"
    grid = np.asarray(grid, dtype=np.uint8)
    if grid.ndim == 4:
        grid = grid[..., 0]
    occupied = int(grid.sum())
    stages.append(Stage("Voxelize", time.perf_counter() - start, f"{backend}, occupied {occupied:,}, fill {occupied / total_cells:.8f}"))

    suffix = f"r{args.resolution}" if args.resolution else f"target{args.target_cells}"
    spacing = tuple(spans[i] / max(1, [nx, ny, nz][i]) for i in range(3))

    occupancy_path = out_dir / f"{input_path.stem}_occupancy_{suffix}_{nx}x{ny}x{nz}.vtk"
    start = time.perf_counter()
    write_occupancy_vtk(occupancy_path, grid, mins, spacing)
    stages.append(Stage("Write Occupancy VTK", time.perf_counter() - start, f"{occupancy_path.name}, {occupancy_path.stat().st_size / 1048576:.2f} MB"))

    if args.surface:
        surface_path = out_dir / f"{input_path.stem}_surface_{suffix}_{nx}x{ny}x{nz}.vtk"
        start = time.perf_counter()
        face_count, point_count = write_surface_vtk(surface_path, grid, mins, spacing)
        stages.append(Stage("Extract Surface VTK", time.perf_counter() - start, f"{face_count:,} quad faces, {point_count:,} points, {surface_path.stat().st_size / 1048576:.2f} MB"))

    stages.append(Stage("Total", time.perf_counter() - total_start, "all stages"))
    return stages


def print_report(stages: list[Stage]) -> None:
    total = stages[-1].seconds if stages and stages[-1].name == "Total" else sum(s.seconds for s in stages)
    cumulative = 0.0
    print("\nTime comparison")
    print("| Stage | Seconds | Percent | Cumulative | Detail |")
    print("|---|---:|---:|---:|---|")
    for stage in stages:
        if stage.name != "Total":
            cumulative += stage.seconds
        pct = stage.seconds / total * 100.0 if total else 0.0
        cum = total if stage.name == "Total" else cumulative
        print(f"| {stage.name} | {stage.seconds:.3f} | {pct:.1f}% | {cum:.3f} | {stage.detail} |")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert STEP/STP/STL to voxel VTK.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "output")
    parser.add_argument("--target-cells", type=parse_positive_int, default=100_000_000)
    parser.add_argument("--resolution", type=parse_positive_int, default=None, help="Override longest-axis voxel count.")
    parser.add_argument("--mesh-size", type=float, default=None)
    parser.add_argument("--ray", default="xyz")
    parser.add_argument("--reuse-stl", action="store_true")
    parser.add_argument("--python", action="store_true", help="Use pure Python voxelizer instead of C++.")
    parser.add_argument("--surface", action="store_true", help="Also write directly visible exterior surface VTK.")
    args = parser.parse_args()

    try:
        stages = run(args)
    except StepMeshingError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2) from None
    print_report(stages)


if __name__ == "__main__":
    main()
