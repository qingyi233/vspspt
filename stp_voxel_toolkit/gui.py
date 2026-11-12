"""Unified GUI for generating voxel models from CAD or occupancy VTK.

Pipeline:
1. STEP/STP -> triangulated STL with Gmsh.
2. STL -> structured occupancy grid with pyvoxel/cvoxelise.
3. Occupancy grid -> legacy binary VTK STRUCTURED_POINTS.
4. Optional occupancy grid -> visible exterior surface VTK UNSTRUCTURED_GRID.
"""

from __future__ import annotations

import math
import queue
import os
import shutil
import sys
import threading
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np

from stl_surface_to_vgm import surface_stl_to_voxels
from vtk_selection_to_vgm import iter_solid_load_voxels, read_occupancy, selected_voxels, write_vgm_stream
from vgm_stats import format_statistics
from voxel_grid import mesh_bounds, resolve_grid, validate_box
from run_dlutfem import DEFAULT_SOLVER_BIN, DEFAULT_WORK_ROOT, run_dlutfem, safe_job_name


ROOT = Path(__file__).resolve().parent
FORCE_PER_CELL = "每细单元力"
FORCE_TOTAL = "选区总力"
PYVOXEL_DIR = ROOT / "pyvoxel_runtime"
if str(PYVOXEL_DIR) not in sys.path:
    sys.path.insert(0, str(PYVOXEL_DIR))

try:
    import gmsh
    import pvoxel
    import voxel
except Exception as exc:  # surfaced in GUI on startup
    gmsh = None
    pvoxel = None
    voxel = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


@dataclass
class RunConfig:
    step_path: Path
    out_dir: Path | None
    target_cells: int
    longest_resolution: int | None
    mesh_size: float | None
    ray: str
    reuse_stl: bool
    use_cpp: bool
    write_occupancy: bool
    write_surface: bool
    box_origin: tuple[float, float, float] | None = None
    box_size: tuple[float, float, float] | None = None
    write_vgm: bool = False
    fine: tuple[int, int, int] = (5, 5, 5)
    vgm_options: dict | None = None


class StepMeshingError(RuntimeError):
    pass


def parse_positive_int(text: str, field_name: str) -> int:
    cleaned = text.strip().replace(",", "")
    if not cleaned:
        raise ValueError(f"{field_name}不能为空")
    value = float(cleaned)
    if not value.is_integer():
        raise ValueError(f"{field_name}必须是整数")
    result = int(value)
    if result <= 0:
        raise ValueError(f"{field_name}必须大于 0")
    return result


def parse_vector3(text: str, field_name: str) -> tuple[float, float, float]:
    values = [part for part in text.replace("，", ",").replace(" ", ",").split(",") if part]
    if len(values) != 3:
        raise ValueError(f"{field_name}必须是 3 个数字，例如 0,0,-1")
    return float(values[0]), float(values[1]), float(values[2])


def log_put(q: queue.Queue, message: str) -> None:
    q.put(("log", message))


def stage_put(q: queue.Queue, stage: str, seconds: float, detail: str) -> None:
    q.put(("stage", stage, seconds, detail))


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


def step_to_stl(step_path: Path, stl_path: Path, mesh_size: float | None, q: queue.Queue) -> int:
    with TemporaryDirectory(prefix="stp_voxel_gui_") as temp_name:
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
                    "STEP 转 STL 失败：Gmsh 无法给这个 CAD 模型生成表面三角网格。\n"
                    f"原始错误：{exc}\n\n"
                    "这通常是 STEP 模型中存在周期面、坏面、退化边或缝合问题，不是体素数或 C++ 体素化的问题。\n"
                    "处理建议：\n"
                    "1. 在 CAD 软件中修复/缝合几何后重新导出 STEP。\n"
                    "2. 直接从 CAD 软件导出 STL，然后在本工具中选择 STL 文件继续体素化。\n"
                    "3. 如果有多个零件，先拆分或简化模型后再转换。"
                ) from exc

            element_types, element_tags, _ = gmsh.model.mesh.getElements(dim=2)
            element_count = sum(tags.size for tags in element_tags)
            log_put(q, f"Gmsh 表面单元分组: {len(element_types)}")
            log_put(q, f"Gmsh 表面单元数: {element_count:,}")
            if element_count == 0:
                raise RuntimeError("Gmsh 没有生成表面单元")

            triangles = collect_gmsh_triangles()
            write_ascii_stl(temp_stl, temp_step.stem, triangles)
        finally:
            gmsh.finalize()

        if not temp_stl.exists() or temp_stl.stat().st_size == 0:
            raise RuntimeError("STL 写出结果为空")
        shutil.copy2(temp_stl, stl_path)
        return len(triangles) // 3


def write_occupancy_vtk(
    path: Path,
    occupancy: np.ndarray,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
) -> None:
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
    masks: list[tuple[str, np.ndarray]] = []

    m = occ.copy()
    m[1:, :, :] &= ~occ[:-1, :, :]
    masks.append(("xmin", m))

    m = occ.copy()
    m[:-1, :, :] &= ~occ[1:, :, :]
    masks.append(("xmax", m))

    m = occ.copy()
    m[:, 1:, :] &= ~occ[:, :-1, :]
    masks.append(("ymin", m))

    m = occ.copy()
    m[:, :-1, :] &= ~occ[:, 1:, :]
    masks.append(("ymax", m))

    m = occ.copy()
    m[:, :, 1:] &= ~occ[:, :, :-1]
    masks.append(("zmin", m))

    m = occ.copy()
    m[:, :, :-1] &= ~occ[:, :, 1:]
    masks.append(("zmax", m))
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
    n = i.size
    pts = np.empty((n, 4, 3), dtype=">f4")

    if direction == "xmin":
        pts[:, 0, :] = np.column_stack([xi, yj, zk])
        pts[:, 1, :] = np.column_stack([xi, yj, zkp])
        pts[:, 2, :] = np.column_stack([xi, yjp, zkp])
        pts[:, 3, :] = np.column_stack([xi, yjp, zk])
    elif direction == "xmax":
        pts[:, 0, :] = np.column_stack([xip, yj, zk])
        pts[:, 1, :] = np.column_stack([xip, yjp, zk])
        pts[:, 2, :] = np.column_stack([xip, yjp, zkp])
        pts[:, 3, :] = np.column_stack([xip, yj, zkp])
    elif direction == "ymin":
        pts[:, 0, :] = np.column_stack([xi, yj, zk])
        pts[:, 1, :] = np.column_stack([xip, yj, zk])
        pts[:, 2, :] = np.column_stack([xip, yj, zkp])
        pts[:, 3, :] = np.column_stack([xi, yj, zkp])
    elif direction == "ymax":
        pts[:, 0, :] = np.column_stack([xi, yjp, zk])
        pts[:, 1, :] = np.column_stack([xi, yjp, zkp])
        pts[:, 2, :] = np.column_stack([xip, yjp, zkp])
        pts[:, 3, :] = np.column_stack([xip, yjp, zk])
    elif direction == "zmin":
        pts[:, 0, :] = np.column_stack([xi, yj, zk])
        pts[:, 1, :] = np.column_stack([xi, yjp, zk])
        pts[:, 2, :] = np.column_stack([xip, yjp, zk])
        pts[:, 3, :] = np.column_stack([xip, yj, zk])
    elif direction == "zmax":
        pts[:, 0, :] = np.column_stack([xi, yj, zkp])
        pts[:, 1, :] = np.column_stack([xip, yj, zkp])
        pts[:, 2, :] = np.column_stack([xip, yjp, zkp])
        pts[:, 3, :] = np.column_stack([xi, yjp, zkp])
    else:
        raise ValueError(direction)
    return pts.reshape(-1, 3)


def write_surface_vtk(path: Path, grid: np.ndarray, origin, spacing, q: queue.Queue) -> tuple[int, int]:
    masks = boundary_masks(grid)
    counts = [(direction, int(mask.sum())) for direction, mask in masks]
    face_count = sum(count for _, count in counts)
    point_count = face_count * 4
    nx, ny, nz = grid.shape
    log_put(q, "表面面片: " + ", ".join(f"{d}={c:,}" for d, c in counts))

    def write_int_scalar(f, name: str, values: np.ndarray) -> None:
        f.write(f"SCALARS {name} int 1\nLOOKUP_TABLE default\n".encode("ascii"))
        f.write(np.asarray(values, dtype=">i4").tobytes())
        f.write(b"\n")

    with path.open("wb") as f:
        f.write(
            (
                "# vtk DataFile Version 3.0\n"
                "Occupied voxel exterior surface\n"
                "BINARY\n"
                "DATASET UNSTRUCTURED_GRID\n"
                f"POINTS {point_count} float\n"
            ).encode("ascii")
        )
        for direction, mask in masks:
            if mask.any():
                f.write(face_points(direction, np.nonzero(mask), origin, spacing).tobytes())
        f.write(b"\n")

        f.write(f"CELLS {face_count} {face_count * 5}\n".encode("ascii"))
        offset = 0
        chunk_faces = 1_000_000
        remaining = face_count
        while remaining:
            n = min(chunk_faces, remaining)
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
            n = min(chunk_faces, remaining)
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


def map_selection(path, occupancy, origin, spacing, label, q, *, workers=1, solid_only=False):
    if path is None:
        return set()
    q.put(("status", f"正在读取{label}..."))
    if path.suffix.lower() == ".stl":
        result = surface_stl_to_voxels(
            path, occupancy, origin, spacing, math.hypot(*spacing), workers=workers,
            progress=lambda done, total, count: log_put(q, f"{label}: {done:,}/{total:,} 三角面，命中 {count:,} 体素"))
    else:
        result = selected_voxels(path, occupancy.shape)
        if solid_only:
            result = set(iter_solid_load_voxels(occupancy, result))
    unit = "个实体细单元" if solid_only else "个体素"
    log_put(q, f"{label}: {path}，命中 {len(result):,} {unit}")
    return result


def write_grid_vgm(occupancy, origin, spacing, output_path, options, q):
    """Use the same boundary mapping and VGM output for CAD and occupancy inputs."""
    workers = options.get("workers", 1)
    constraints = map_selection(options.get("constraint"), occupancy, origin, spacing, "约束选区", q, workers=workers)
    loads = map_selection(options.get("load"), occupancy, origin, spacing, "载荷选区", q,
                          workers=workers, solid_only=True)
    force = options.get("force", 0.)
    if options.get("force_mode") == FORCE_TOTAL:
        if not loads:
            raise ValueError("载荷选区未命中实体细单元，无法分配总力")
        total_force = force
        force /= len(loads)
        log_put(q, f"选区总力 {total_force:.8g} 均分至 {len(loads):,} 个实体细单元，每单元 {force:.8g}")
    q.put(("status", "正在写出 VGM 并统计网格..."))
    stats = write_vgm_stream(
        output_path, occupancy, origin, spacing, constraints, loads, options["fine"],
        options.get("spc_dof", "123"), 1, force,
        options.get("direction", (0., 0., -1.)), 1, "Voxel model",
        progress=lambda stage, done, total: log_put(q, f"{stage}: {done:,}/{total:,}"),
        spc_node_mode=options.get("spc_node_mode", "cell-corners"))
    log_put(q, f"VGM 输出: {output_path}")
    log_put(q, f"统计文件: {output_path.with_suffix('.stats.json')}")
    q.put(("vgm_statistics", format_statistics(stats)))
    q.put(("solve_vgm_path", str(output_path)))


def voxelize_model(config: RunConfig, q: queue.Queue, stl_dir: Path):
    """Build the same occupancy grid for preview and file generation."""
    input_path = config.step_path
    stl_name = f"{input_path.stem}.stl"
    stl_path = stl_dir / stl_name
    cached_stl = config.out_dir / stl_name if config.out_dir is not None else None

    start = time.time()
    if input_path.suffix.lower() == ".stl":
        stl_path = input_path
        stage_put(q, "STEP 转 STL", 0.0, "输入已是 STL，跳过")
    elif input_path.suffix.lower() in {".stp", ".step"}:
        if config.reuse_stl and cached_stl and cached_stl.is_file() and cached_stl.stat().st_size > 0:
            stl_path = cached_stl
            stage_put(q, "STEP 转 STL", 0.0, f"复用 {stl_path}")
        else:
            triangles = step_to_stl(input_path, stl_path, config.mesh_size, q)
            stage_put(q, "STEP 转 STL", time.time() - start, f"{triangles:,} 个三角面, {stl_path.stat().st_size / 1048576:.2f} MB")
    else:
        raise ValueError("输入文件必须是 .stp、.step 或 .stl")

    start = time.time()
    mesh_vertices = pvoxel.Read_stl(str(stl_path))[0]
    if mesh_vertices is None or mesh_vertices.size == 0:
        raise RuntimeError(f"没有从 {stl_path} 读取到三角面")
    (nx, ny, nz), mins, spacing = resolve_grid(
        mesh_vertices, config.target_cells, config.longest_resolution, config.box_origin, config.box_size,
        fine=config.fine if config.write_vgm else None)
    maxs = mesh_bounds(mesh_vertices)[1] if config.box_size is None else np.asarray(mins) + config.box_size
    bounds = np.array([mins, maxs])
    total_cells = nx * ny * nz
    detail = f"{mesh_vertices.shape[0]:,} 个三角面, 网格 {nx}x{ny}x{nz}={total_cells:,}"
    stage_put(q, "读取 STL / 生成网格", time.time() - start, detail)
    log_put(q, f"包围盒最小值: {mins}")
    log_put(q, f"包围盒最大值: {maxs}")
    log_put(q, f"预计占据场内存: {total_cells / 1e6:.1f} MB")

    start = time.time()
    if config.use_cpp:
        vox_obj = voxel.Voxel(nx, ny, nz)
        grid = vox_obj.VOXELISE(mesh_vertices, ray=config.ray, parallel=False, bounds=bounds)
    else:
        vox_obj = pvoxel.Voxel(nx, ny, nz)
        grid = vox_obj.VOXELISE(mesh_vertices, ray=config.ray, parallel=False, bounds=bounds)
    grid = np.asarray(grid, dtype=np.uint8)
    if grid.ndim == 4:
        grid = grid[..., 0]
    occupied = int(np.count_nonzero(grid))
    stage_put(q, "体素化", time.time() - start, f"占据体素 {occupied:,}, 填充率 {occupied / total_cells:.8f}")
    return grid, mins, spacing, occupied


def run_pipeline(config: RunConfig, q: queue.Queue) -> None:
    total_start = time.time()
    if config.out_dir is None:
        raise ValueError("请选择输出目录")
    config.out_dir.mkdir(parents=True, exist_ok=True)
    input_path = config.step_path
    suffix = f"r{config.longest_resolution}" if config.longest_resolution else f"target{config.target_cells}"
    log_put(q, f"输入文件: {input_path}")
    log_put(q, f"输出目录: {config.out_dir}")
    grid, mins, spacing, _ = voxelize_model(config, q, config.out_dir)
    nx, ny, nz = grid.shape

    size_tag = f"{nx}x{ny}x{nz}"
    occupancy_path = config.out_dir / f"{input_path.stem}_occupancy_{suffix}_{size_tag}.vtk"
    surface_path = config.out_dir / f"{input_path.stem}_surface_{suffix}_{size_tag}.vtk"

    if config.write_occupancy:
        start = time.time()
        write_occupancy_vtk(occupancy_path, grid, mins, spacing)
        stage_put(q, "写出占据场 VTK", time.time() - start, f"{occupancy_path.name}, {occupancy_path.stat().st_size / 1048576:.2f} MB")
        log_put(q, f"占据场 VTK: {occupancy_path}")

    if config.write_surface:
        start = time.time()
        face_count, point_count = write_surface_vtk(surface_path, grid, mins, spacing, q)
        detail = f"{face_count:,} 个四边形面, {point_count:,} 个点, {surface_path.stat().st_size / 1048576:.2f} MB"
        stage_put(q, "提取表面 VTK", time.time() - start, detail)
        log_put(q, f"可视表面 VTK: {surface_path}")

    if config.write_vgm:
        start = time.time()
        vgm_path = config.out_dir / f"{input_path.stem}.vgm"
        write_grid_vgm(grid, mins, spacing, vgm_path,
                       config.vgm_options or {"fine": config.fine}, q)
        stage_put(q, "写出 VGM / 统计", time.time() - start, vgm_path.name)

    stage_put(q, "总计", time.time() - total_start, "完成")
    q.put(("done", config.write_vgm))


class VoxelGui(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("体素模型工具")
        self.geometry("1020x800")
        self.minsize(880, 640)
        self.queue: queue.Queue = queue.Queue()
        self.worker: threading.Thread | None = None
        self.step_var = tk.StringVar()
        self.out_var = tk.StringVar(value=str(ROOT / "output" / "gui_voxel_output"))
        self.grid_mode_var = tk.StringVar(value="目标体素数")
        self.target_var = tk.StringVar(value="100000000")
        self.resolution_var = tk.StringVar(value="400")
        self.mesh_size_var = tk.StringVar()
        self.ray_var = tk.StringVar(value="xyz")
        self.reuse_stl_var = tk.BooleanVar(value=True)
        self.use_cpp_var = tk.BooleanVar(value=True)
        self.write_occ_var = tk.BooleanVar(value=False)
        self.write_surface_var = tk.BooleanVar(value=False)
        self.write_vgm_var = tk.BooleanVar(value=True)
        self.custom_box_var = tk.BooleanVar(value=False)
        self.box_origin_var = tk.StringVar(value="0,0,0")
        self.box_size_var = tk.StringVar(value="100,100,100")
        self.vgm_stats_var = tk.StringVar()
        self.vgm_constraint_var = tk.StringVar()
        self.vgm_load_var = tk.StringVar()
        self.vgm_fine_x_var = tk.StringVar(value="5")
        self.vgm_fine_y_var = tk.StringVar(value="5")
        self.vgm_fine_z_var = tk.StringVar(value="5")
        self.vgm_fine_preset_var = tk.StringVar(value="5x5x5")
        self.vgm_spc_dof_var = tk.StringVar(value="123")
        self.vgm_force_mode_var = tk.StringVar(value=FORCE_PER_CELL)
        self.vgm_force_var = tk.StringVar(value="1000")
        self.vgm_direction_var = tk.StringVar(value="0,0,-1")
        self.vgm_workers_var = tk.StringVar(value=str(max(1, min(8, (os.cpu_count() or 2) // 2))))
        self.vgm_spc_node_mode_var = tk.StringVar(value="命中粗单元8角点")
        self.solve_vgm_var = tk.StringVar()
        self.solve_out_var = tk.StringVar(value=str(ROOT / "output" / "dlutfem_result"))
        self.solve_job_var = tk.StringVar()
        self.solve_solver_bin_var = tk.StringVar(value=str(DEFAULT_SOLVER_BIN))
        self.solve_work_root_var = tk.StringVar(value=str(DEFAULT_WORK_ROOT))
        self.solve_refresh_bin_var = tk.BooleanVar(value=False)
        self.grid_count_var = tk.StringVar()
        self.load_count_var = tk.StringVar()
        self.grid_count_revision = 0
        self.load_count_revision = 0
        self.status_var = tk.StringVar(value="就绪")
        self._build_ui()
        for var in (self.vgm_fine_x_var, self.vgm_fine_y_var, self.vgm_fine_z_var):
            var.trace_add("write", self._sync_fine_preset)
        for var in (self.step_var, self.write_vgm_var, self.custom_box_var,
                    self.vgm_constraint_var, self.vgm_load_var, self.vgm_stats_var, self.vgm_fine_preset_var):
            var.trace_add("write", self._refresh_generation_ui)
        for var in (self.step_var, self.out_var, self.grid_mode_var, self.target_var, self.resolution_var,
                    self.mesh_size_var, self.ray_var, self.reuse_stl_var, self.use_cpp_var, self.write_vgm_var,
                    self.custom_box_var, self.box_origin_var, self.box_size_var,
                    self.vgm_fine_x_var, self.vgm_fine_y_var, self.vgm_fine_z_var):
            var.trace_add("write", self._invalidate_grid_count)
        self.vgm_load_var.trace_add("write", self._invalidate_load_count)
        self._refresh_generation_ui()
        self.after(100, self._poll_queue)

    def _scroll_page(self, notebook, title):
        page = ttk.Frame(notebook)
        notebook.add(page, text=title)
        canvas = tk.Canvas(page, highlightthickness=0, height=440)
        scrollbar = ttk.Scrollbar(page, orient="vertical", command=canvas.yview)
        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        canvas.configure(yscrollcommand=scrollbar.set)
        body = ttk.Frame(canvas, padding=(0, 8, 8, 0))
        window = canvas.create_window((0, 0), window=body, anchor="nw")
        body.bind("<Configure>", lambda _event: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda event: canvas.itemconfigure(window, width=event.width))

        def scroll(event):
            widget = self.winfo_containing(event.x_root, event.y_root)
            while widget is not None:
                if widget == page:
                    if body.winfo_reqheight() > canvas.winfo_height():
                        canvas.yview_scroll(int(-event.delta / 120), "units")
                    return
                widget = widget.master
        self.bind("<MouseWheel>", scroll, add="+")
        return body

    def _advanced_group(self, parent):
        wrapper = ttk.Frame(parent)
        wrapper.pack(fill="x", pady=(8, 0))
        opened = tk.BooleanVar(value=False)
        body = ttk.Frame(wrapper, padding=(8, 6, 0, 0))
        ttk.Checkbutton(wrapper, text="高级设置", variable=opened,
                        command=lambda: self._show(body, opened.get())).pack(anchor="w")
        return wrapper, body

    @staticmethod
    def _show(widget, visible):
        if visible:
            widget.pack(fill="x", pady=(6, 0))
        else:
            widget.pack_forget()

    def _path_row(self, parent, row, label, variable, browse):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=4)
        entry = ttk.Entry(parent, textvariable=variable)
        entry.grid(row=row, column=1, sticky="ew", padx=8, pady=4)
        ttk.Button(parent, text="浏览", command=browse).grid(row=row, column=2, pady=4)
        parent.columnconfigure(1, weight=1)
        return entry

    def _build_ui(self) -> None:
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)
        main.columnconfigure(0, weight=1)
        main.rowconfigure(0, weight=1)
        self.notebook = ttk.Notebook(main)
        self.notebook.grid(row=0, column=0, sticky="nsew")
        generation = self._scroll_page(self.notebook, "模型生成")
        self.generation_canvas = generation.master
        solve = self._scroll_page(self.notebook, "有限元求解")

        files = ttk.LabelFrame(generation, text="输入与输出", padding=10)
        files.pack(fill="x")
        self._path_row(files, 0, "模型 / 占据场", self.step_var, self._browse_step)
        self._path_row(files, 1, "输出目录", self.out_var, lambda: self._browse_directory(self.out_var))
        ttk.Label(files, text="支持 STEP / STL / occupancy VTK；VGM 按输入文件名自动命名。",
                  foreground="#555555").grid(row=2, column=0, columnspan=3, sticky="w", pady=(4, 0))

        self.grid_frame = ttk.LabelFrame(generation, text="体素网格", padding=10)
        self.grid_frame.pack(fill="x", pady=(8, 0))
        row = ttk.Frame(self.grid_frame)
        row.pack(fill="x")
        ttk.Label(row, text="分辨率方式").pack(side="left")
        mode = ttk.Combobox(row, textvariable=self.grid_mode_var,
                            values=("目标体素数", "最长轴格数"), state="readonly", width=14)
        mode.pack(side="left", padx=8)
        mode.bind("<<ComboboxSelected>>", self._grid_mode_changed)
        self.grid_value_entry = ttk.Entry(row, textvariable=self.target_var, width=18)
        self.grid_value_entry.pack(side="left")
        ttk.Checkbutton(row, text="自定义包围盒", variable=self.custom_box_var).pack(side="left", padx=(24, 0))
        self.box_fields = ttk.Frame(self.grid_frame)
        ttk.Label(self.box_fields, text="最小角点 X,Y,Z").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.box_fields, textvariable=self.box_origin_var, width=24).grid(row=0, column=1, padx=8)
        ttk.Label(self.box_fields, text="尺寸 X,Y,Z").grid(row=0, column=2, padx=(12, 0))
        ttk.Entry(self.box_fields, textvariable=self.box_size_var, width=24).grid(row=0, column=3, padx=8)
        ttk.Label(self.box_fields, text="单位同模型；未完整包住模型时会报错停止。", foreground="#555555").grid(
            row=1, column=0, columnspan=4, sticky="w", pady=(4, 0))

        self.outputs_frame = ttk.Frame(generation)
        self.outputs_frame.pack(fill="x", pady=(8, 0))
        ttk.Label(self.outputs_frame, text="输出").pack(side="left", padx=(0, 8))
        ttk.Checkbutton(self.outputs_frame, text="VGM", variable=self.write_vgm_var).pack(side="left", padx=(0, 16))
        ttk.Checkbutton(self.outputs_frame, text="占据场 VTK", variable=self.write_occ_var).pack(side="left", padx=(0, 16))
        ttk.Checkbutton(self.outputs_frame, text="表面 VTK", variable=self.write_surface_var).pack(side="left")

        self.vgm_frame = ttk.LabelFrame(generation, text="VGM 设置", padding=10)
        self.vgm_frame.pack(fill="x", pady=(8, 0))
        row = ttk.Frame(self.vgm_frame)
        row.pack(fill="x")
        ttk.Label(row, text="粗单元细分").pack(side="left")
        preset = ttk.Combobox(row, textvariable=self.vgm_fine_preset_var,
                             values=("5x5x5", "10x10x10", "自定义"), state="readonly", width=12)
        preset.pack(side="left", padx=8)
        preset.bind("<<ComboboxSelected>>", self._apply_fine_preset)
        self.fine_fields = ttk.Frame(row)
        for label, var in zip(("X", "Y", "Z"),
                              (self.vgm_fine_x_var, self.vgm_fine_y_var, self.vgm_fine_z_var)):
            ttk.Label(self.fine_fields, text=label).pack(side="left")
            ttk.Entry(self.fine_fields, textvariable=var, width=5).pack(side="left", padx=(4, 8))
        self.fine_hint = ttk.Label(self.vgm_frame, foreground="#555555")
        selections = ttk.Frame(self.vgm_frame)
        selections.pack(fill="x", pady=(4, 0))
        self._path_row(selections, 0, "约束选区（可选）", self.vgm_constraint_var,
                       lambda: self._browse_selection(self.vgm_constraint_var))
        self._path_row(selections, 1, "载荷选区（可选）", self.vgm_load_var,
                       lambda: self._browse_selection(self.vgm_load_var))
        self.constraint_fields = ttk.Frame(self.vgm_frame)
        ttk.Label(self.constraint_fields, text="约束自由度").pack(side="left")
        ttk.Entry(self.constraint_fields, textvariable=self.vgm_spc_dof_var, width=8).pack(side="left", padx=8)
        self.load_fields = ttk.Frame(self.vgm_frame)
        load_preview = ttk.Frame(self.load_fields)
        load_preview.pack(fill="x")
        self.load_count_button = ttk.Button(load_preview, text="计算命中数",
                                            command=lambda: self._start_grid_count(include_load=True))
        self.load_count_button.pack(side="left")
        ttk.Label(load_preview, textvariable=self.load_count_var).pack(side="left", padx=8)
        force_fields = ttk.Frame(self.load_fields)
        force_fields.pack(fill="x", pady=(6, 0))
        ttk.Combobox(force_fields, textvariable=self.vgm_force_mode_var, state="readonly",
                     values=(FORCE_PER_CELL, FORCE_TOTAL), width=12).pack(side="left")
        ttk.Entry(force_fields, textvariable=self.vgm_force_var, width=12).pack(side="left", padx=8)
        ttk.Label(force_fields, text="方向 X,Y,Z").pack(side="left", padx=(12, 0))
        ttk.Entry(force_fields, textvariable=self.vgm_direction_var, width=18).pack(side="left", padx=8)

        self.generation_advanced, advanced = self._advanced_group(generation)
        self.cad_advanced_fields = ttk.Frame(advanced)
        ttk.Label(self.cad_advanced_fields, text="曲面网格尺寸").grid(row=0, column=0, sticky="w")
        ttk.Entry(self.cad_advanced_fields, textvariable=self.mesh_size_var, width=14).grid(row=0, column=1, padx=8)
        ttk.Label(self.cad_advanced_fields, text="射线方向").grid(row=0, column=2, padx=(16, 0))
        ttk.Combobox(self.cad_advanced_fields, textvariable=self.ray_var,
                     values=("xyz", "x", "y", "z", "xy", "xz", "yz"), width=8, state="readonly").grid(row=0, column=3, padx=8)
        ttk.Checkbutton(self.cad_advanced_fields, text="使用 C++ 体素化", variable=self.use_cpp_var).grid(
            row=1, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Checkbutton(self.cad_advanced_fields, text="复用已有 STL", variable=self.reuse_stl_var).grid(
            row=1, column=2, columnspan=2, sticky="w", pady=(6, 0))
        self.constraint_advanced_fields = ttk.Frame(advanced)
        ttk.Label(self.constraint_advanced_fields, text="粗节点约束映射").pack(side="left")
        ttk.Combobox(self.constraint_advanced_fields, textvariable=self.vgm_spc_node_mode_var,
                     values=("命中粗单元8角点", "命中粗单元最小角点"), width=22, state="readonly").pack(side="left", padx=8)
        self.stl_advanced_fields = ttk.Frame(advanced)
        ttk.Label(self.stl_advanced_fields, text="曲面选区映射线程").pack(side="left")
        ttk.Entry(self.stl_advanced_fields, textvariable=self.vgm_workers_var, width=8).pack(side="left", padx=8)

        self.stats_frame = ttk.LabelFrame(generation, text="网格统计", padding=10)
        ttk.Label(self.stats_frame, textvariable=self.vgm_stats_var, wraplength=820).pack(anchor="w")

        solver_files = ttk.LabelFrame(solve, text="求解文件", padding=10)
        solver_files.pack(fill="x")
        self._path_row(solver_files, 0, "输入 VGM", self.solve_vgm_var, self._browse_solve_vgm)
        self._path_row(solver_files, 1, "结果目录", self.solve_out_var,
                       lambda: self._browse_directory(self.solve_out_var))
        _, advanced = self._advanced_group(solve)
        self._path_row(advanced, 0, "求解器目录", self.solve_solver_bin_var,
                       lambda: self._browse_directory(self.solve_solver_bin_var))
        self._path_row(advanced, 1, "临时工作目录", self.solve_work_root_var,
                       lambda: self._browse_directory(self.solve_work_root_var))
        ttk.Label(advanced, text="任务名称").grid(row=2, column=0, sticky="w", pady=4)
        ttk.Entry(advanced, textvariable=self.solve_job_var, width=28).grid(row=2, column=1, sticky="w", padx=8)
        ttk.Checkbutton(advanced, text="刷新求解器运行库", variable=self.solve_refresh_bin_var).grid(
            row=3, column=0, columnspan=3, sticky="w", pady=(6, 0))

        actions = ttk.Frame(main)
        actions.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        self.generation_actions = ttk.Frame(actions)
        self.generation_actions.pack(side="left")
        self.run_button = ttk.Button(self.generation_actions, text="开始生成", command=self._start)
        self.run_button.pack(side="left")
        self.grid_count_button = ttk.Button(self.generation_actions, text="统计实体体素", command=self._start_grid_count)
        self.grid_count_button.pack(side="left", padx=(8, 0))
        ttk.Label(self.generation_actions, textvariable=self.grid_count_var).pack(side="left", padx=8)
        self.solve_button = ttk.Button(actions, text="开始求解", command=self._start_solve)
        ttk.Label(actions, textvariable=self.status_var).pack(side="right")
        self.notebook.bind("<<NotebookTabChanged>>", self._tab_changed)
        log_frame = ttk.LabelFrame(main, text="运行日志", padding=8)
        log_frame.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        self.log_text = tk.Text(log_frame, height=7, wrap="word")
        self.log_text.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=scrollbar.set)

    def _tab_changed(self, _event=None):
        generation = self.notebook.index("current") == 0
        self.generation_actions.pack_forget()
        self.solve_button.pack_forget()
        (self.generation_actions if generation else self.solve_button).pack(side="left")

    def _invalidate_grid_count(self, variable, *_args):
        if variable != str(self.step_var) and Path(self.step_var.get().strip()).suffix.lower() == ".vtk":
            return
        self.grid_count_revision += 1
        self.grid_count_var.set("")
        self._invalidate_load_count()

    def _invalidate_load_count(self, *_args):
        self.load_count_revision += 1
        self.load_count_var.set("")

    def _grid_mode_changed(self, _event=None):
        var = self.target_var if self.grid_mode_var.get() == "目标体素数" else self.resolution_var
        self.grid_value_entry.configure(textvariable=var)

    def _refresh_generation_ui(self, *_args):
        is_vtk = Path(self.step_var.get().strip()).suffix.lower() == ".vtk"
        vgm = is_vtk or self.write_vgm_var.get()
        self._show(self.box_fields, self.custom_box_var.get())
        self._show(self.constraint_fields, bool(self.vgm_constraint_var.get().strip()))
        self._show(self.load_fields, bool(self.vgm_load_var.get().strip()))
        custom = self.vgm_fine_preset_var.get() == "自定义"
        self._show(self.fine_fields, custom)
        # Fine fields sit on the same row as the preset selector.
        if custom:
            self.fine_fields.pack_configure(side="left", fill="none", pady=0)
        self.fine_hint.configure(text="10x10x10 需要支持该细分的求解器和 EN_10 权重。"
                                if self.vgm_fine_preset_var.get() == "10x10x10" else "")
        if self.vgm_fine_preset_var.get() == "10x10x10":
            self.fine_hint.pack(after=self.fine_fields.master, anchor="w", pady=(4, 0))
        else:
            self.fine_hint.pack_forget()
        constraint = vgm and bool(self.vgm_constraint_var.get().strip())
        stl = vgm and any(Path(var.get()).suffix.lower() == ".stl" for var in
                         (self.vgm_constraint_var, self.vgm_load_var))
        self._show(self.cad_advanced_fields, not is_vtk)
        self._show(self.constraint_advanced_fields, constraint)
        self._show(self.stl_advanced_fields, stl)
        # Reinsert conditional sections in a stable order when the input type changes.
        sections = ((self.grid_frame, not is_vtk), (self.outputs_frame, not is_vtk),
                    (self.vgm_frame, vgm), (self.generation_advanced, not is_vtk or constraint or stl),
                    (self.stats_frame, bool(self.vgm_stats_var.get())))
        for widget, _visible in sections:
            widget.pack_forget()
        for widget, visible in sections:
            self._show(widget, visible)

    def _browse_step(self):
        path = filedialog.askopenfilename(filetypes=[("模型或占据场", "*.stp *.step *.stl *.vtk"), ("所有文件", "*.*")])
        if path:
            self.step_var.set(path)

    def _browse_directory(self, variable):
        path = filedialog.askdirectory()
        if path:
            variable.set(path)

    def _browse_selection(self, variable):
        path = filedialog.askopenfilename(filetypes=[("选区 VTK / 曲面 STL", "*.vtk *.stl"), ("所有文件", "*.*")])
        if path:
            variable.set(path)

    def _browse_solve_vgm(self):
        path = filedialog.askopenfilename(filetypes=[("VGM 文件", "*.vgm"), ("所有文件", "*.*")])
        if path:
            self.solve_vgm_var.set(path)

    def _sync_fine_preset(self, *_args):
        label = "x".join(var.get().strip() for var in
                         (self.vgm_fine_x_var, self.vgm_fine_y_var, self.vgm_fine_z_var))
        self.vgm_fine_preset_var.set(label if label in ("5x5x5", "10x10x10") else "自定义")

    def _apply_fine_preset(self, _event=None):
        label = self.vgm_fine_preset_var.get()
        if label in ("5x5x5", "10x10x10"):
            for var, value in zip((self.vgm_fine_x_var, self.vgm_fine_y_var, self.vgm_fine_z_var), label.split("x")):
                var.set(value)

    @staticmethod
    def _selection_path(var, label):
        text = var.get().strip()
        if not text:
            return None
        path = Path(text)
        if not path.is_file() or path.suffix.lower() not in {".vtk", ".stl"}:
            raise ValueError(f"{label}必须是存在的 VTK 或 STL 文件")
        return path

    def _parse_fine(self):
        return tuple(parse_positive_int(var.get(), "细分 " + label) for label, var in zip(
            ("X", "Y", "Z"), (self.vgm_fine_x_var, self.vgm_fine_y_var, self.vgm_fine_z_var)))

    def _parse_boundary_config(self):
        constraint = self._selection_path(self.vgm_constraint_var, "约束选区")
        load = self._selection_path(self.vgm_load_var, "载荷选区")
        fine = self._parse_fine()
        spc_dof = self.vgm_spc_dof_var.get().strip() if constraint else "123"
        if constraint and (not spc_dof or not set(spc_dof).issubset("123")):
            raise ValueError("约束自由度只能包含 1、2、3")
        mode = {"命中粗单元8角点": "cell-corners", "命中粗单元最小角点": "cell-min-corner"}.get(
            self.vgm_spc_node_mode_var.get(), "cell-corners")
        workers = parse_positive_int(self.vgm_workers_var.get(), "曲面映射线程") if any(
            path and path.suffix.lower() == ".stl" for path in (constraint, load)) else 1
        force = float(self.vgm_force_var.get()) if load else 0.
        direction = parse_vector3(self.vgm_direction_var.get(), "力方向") if load else (0., 0., -1.)
        if not np.isfinite(force) or not np.isfinite(direction).all():
            raise ValueError("力大小和方向必须是有限数值")
        force_mode = self.vgm_force_mode_var.get() if load else FORCE_PER_CELL
        if force_mode == FORCE_TOTAL:
            scale = max(map(abs, direction))
            if scale == 0:
                raise ValueError("设定选区总力时，力方向必须为非零向量")
            scaled = tuple(component / scale for component in direction)
            length = math.hypot(*scaled)
            direction = tuple(component / length for component in scaled)
        return dict(constraint=constraint, load=load, fine=fine, spc_dof=spc_dof,
                    spc_node_mode=mode, workers=workers, force=force, direction=direction, force_mode=force_mode)

    def _input_path(self, suffixes):
        path = Path(self.step_var.get().strip())
        if not path.is_file() or path.suffix.lower() not in suffixes:
            raise ValueError("请选择有效的 STEP、STL 或 occupancy VTK 文件")
        if path.suffix.lower() == ".vtk" and "_surface_" in path.name.lower():
            raise ValueError("请选择 occupancy VTK，不能使用 surface VTK")
        return path

    def _output_directory(self):
        text = self.out_var.get().strip()
        if not text:
            raise ValueError("请选择输出目录")
        return Path(text)

    def _parse_config(self):
        step = self._input_path({".stp", ".step", ".stl"})
        out_dir = self._output_directory()
        if not (self.write_occ_var.get() or self.write_surface_var.get() or self.write_vgm_var.get()):
            raise ValueError("至少选择一种输出格式")
        options = self._parse_boundary_config() if self.write_vgm_var.get() else None
        return self._parse_cad_grid(step, out_dir, options)

    def _parse_cad_grid(self, step, out_dir, options):
        resolution = None
        target = 100_000_000
        if self.grid_mode_var.get() == "最长轴格数":
            resolution = parse_positive_int(self.resolution_var.get(), "最长轴格数")
            if resolution < 2:
                raise ValueError("最长轴格数必须大于 1")
        else:
            target = parse_positive_int(self.target_var.get(), "目标体素数")
        mesh_text = self.mesh_size_var.get().strip() if step.suffix.lower() != ".stl" else ""
        mesh_size = float(mesh_text) if mesh_text else None
        if mesh_size is not None and (not np.isfinite(mesh_size) or mesh_size <= 0):
            raise ValueError("曲面网格尺寸必须为有限正数")
        origin = size = None
        if self.custom_box_var.get():
            origin, size = validate_box(parse_vector3(self.box_origin_var.get(), "包围盒最小角点"),
                                        parse_vector3(self.box_size_var.get(), "包围盒尺寸"))
        return RunConfig(step_path=step, out_dir=out_dir, target_cells=target,
                         longest_resolution=resolution, mesh_size=mesh_size, ray=self.ray_var.get(),
                         reuse_stl=self.reuse_stl_var.get(), use_cpp=self.use_cpp_var.get(),
                         write_occupancy=self.write_occ_var.get(), write_surface=self.write_surface_var.get(),
                         box_origin=origin, box_size=size, write_vgm=options is not None,
                         fine=options["fine"] if options else (5, 5, 5), vgm_options=options)

    def _parse_vgm_config(self):
        occupancy = self._input_path({".vtk"})
        return dict(occupancy=occupancy, output=self._output_directory() / f"{occupancy.stem}.vgm",
                    **self._parse_boundary_config())

    def _begin_task(self, target, config, status, clear_statistics=True):
        self.log_text.delete("1.0", "end")
        if clear_statistics:
            self.vgm_stats_var.set("")
        self._set_task_buttons("disabled")
        self.status_var.set(status)
        self.worker = threading.Thread(target=target, args=(config,), daemon=True)
        self.worker.start()

    def _start(self):
        if self.worker and self.worker.is_alive():
            return
        try:
            if Path(self.step_var.get().strip()).suffix.lower() == ".vtk":
                config = self._parse_vgm_config()
                target = self._run_vgm_worker
            else:
                if IMPORT_ERROR is not None:
                    raise ValueError(f"CAD 体素化依赖导入失败: {IMPORT_ERROR}")
                config = self._parse_config()
                target = self._run_worker
        except Exception as exc:
            messagebox.showerror("参数错误", str(exc))
            return
        self._begin_task(target, config, "正在生成模型...")

    def _start_grid_count(self, include_load=False):
        if self.worker and self.worker.is_alive():
            return
        try:
            model = self._input_path({".stp", ".step", ".stl", ".vtk"})
            load = self._selection_path(self.vgm_load_var, "载荷选区") if include_load else None
            if include_load and load is None:
                raise ValueError("请先选择载荷选区")
            workers = 1
            if load is not None and load.suffix.lower() == ".stl":
                workers = parse_positive_int(self.vgm_workers_var.get(), "曲面映射线程")
            if model.suffix.lower() != ".vtk":
                if IMPORT_ERROR is not None:
                    raise ValueError(f"CAD 体素化依赖导入失败: {IMPORT_ERROR}")
                out_text = self.out_var.get().strip()
                options = {"fine": self._parse_fine()} if self.write_vgm_var.get() else None
                model = self._parse_cad_grid(model, Path(out_text) if out_text else None, options)
        except (ValueError, OSError) as exc:
            messagebox.showerror("参数错误", str(exc))
            return
        self.grid_count_var.set("")
        if include_load:
            self.load_count_var.set("")
        request = (model, load, workers, self.grid_count_revision, self.load_count_revision)
        self._begin_task(self._run_grid_count_worker, request, "正在统计体素...", clear_statistics=False)

    def _run_grid_count_worker(self, request):
        model, load, workers, grid_revision, load_revision = request
        try:
            if isinstance(model, RunConfig):
                with TemporaryDirectory(prefix="voxel_count_") as scratch:
                    occupancy, origin, spacing, solid_count = voxelize_model(model, self.queue, Path(scratch))
            else:
                log_put(self.queue, f"读取占据场: {model}")
                occupancy, origin, spacing = read_occupancy(model)
                solid_count = int(np.count_nonzero(occupancy))
            hit_count = None
            if load is not None:
                hits = map_selection(load, occupancy, origin, spacing, "载荷选区", self.queue,
                                     workers=workers, solid_only=True)
                hit_count = len(hits)
            self.queue.put(("grid_counts", grid_revision, load_revision, solid_count, hit_count))
        except Exception as exc:
            self.queue.put(("error", str(exc), traceback.format_exc()))

    def _parse_solve_config(self):
        vgm_path = Path(self.solve_vgm_var.get().strip())
        solver_bin = Path(self.solve_solver_bin_var.get().strip())
        work_root = Path(self.solve_work_root_var.get().strip())
        if not vgm_path.is_file() or vgm_path.suffix.lower() != ".vgm":
            raise ValueError("请选择有效的 VGM 文件")
        if not self.solve_out_var.get().strip():
            raise ValueError("请选择结果目录")
        if not (solver_bin / "DLUTFEM.exe").is_file():
            raise ValueError("求解器目录中找不到 DLUTFEM.exe")
        if not self.solve_work_root_var.get().strip() or not str(work_root).isascii():
            raise ValueError("临时工作目录不能为空或包含非 ASCII 字符")
        return dict(vgm=vgm_path, output_dir=Path(self.solve_out_var.get().strip()), solver_bin=solver_bin,
                    work_root=work_root, job=self.solve_job_var.get().strip() or safe_job_name(vgm_path.stem),
                    refresh_bin=self.solve_refresh_bin_var.get())

    def _start_solve(self):
        if self.worker and self.worker.is_alive():
            return
        try:
            config = self._parse_solve_config()
        except Exception as exc:
            messagebox.showerror("参数错误", str(exc))
            return
        self._begin_task(self._run_solve_worker, config, "正在求解...", clear_statistics=False)

    def _run_worker(self, config):
        try:
            run_pipeline(config, self.queue)
        except Exception as exc:
            self.queue.put(("error", str(exc), traceback.format_exc()))

    def _run_vgm_worker(self, config):
        try:
            start = time.time()
            log_put(self.queue, f"读取占据场: {config['occupancy']}")
            occupancy, origin, spacing = read_occupancy(config["occupancy"])
            write_grid_vgm(occupancy, origin, spacing, config["output"], config, self.queue)
            stage_put(self.queue, "总计", time.time() - start, "VGM 已生成")
            self.queue.put(("done", True))
        except Exception as exc:
            self.queue.put(("error", str(exc), traceback.format_exc()))

    def _run_solve_worker(self, config):
        try:
            start = time.time()
            returncode = run_dlutfem(config["vgm"], config["solver_bin"], config["work_root"],
                                     config["output_dir"], config["job"], bool(config["refresh_bin"]),
                                     log=lambda message: log_put(self.queue, message))
            self.queue.put(("solve_done", time.time() - start, returncode, str(config["output_dir"])))
        except Exception as exc:
            self.queue.put(("error", str(exc), traceback.format_exc()))

    def _set_task_buttons(self, state):
        self.run_button.configure(state=state)
        self.solve_button.configure(state=state)
        self.grid_count_button.configure(state=state)
        self.load_count_button.configure(state=state)

    def _append_log(self, message):
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")

    def _poll_queue(self):
        try:
            while True:
                item = self.queue.get_nowait()
                kind = item[0]
                if kind == "log":
                    self._append_log(item[1])
                elif kind == "status":
                    self.status_var.set(item[1])
                elif kind == "stage":
                    _, name, seconds, detail = item
                    self._append_log(f"{name}: {seconds:.3f} 秒；{detail}")
                elif kind == "vgm_statistics":
                    self.vgm_stats_var.set(item[1])
                    self.update_idletasks()
                    self.generation_canvas.yview_moveto(1.0)
                elif kind == "solve_vgm_path":
                    self.solve_vgm_var.set(item[1])
                elif kind == "grid_counts":
                    _, grid_revision, load_revision, solid_count, hit_count = item
                    if grid_revision != self.grid_count_revision:
                        self.status_var.set("网格设置已变更，请重新统计")
                    else:
                        self.grid_count_var.set(f"实体体素 {solid_count:,}")
                        if hit_count is not None and load_revision != self.load_count_revision:
                            self.status_var.set("载荷选区已变更，请重新计算命中数")
                        else:
                            if hit_count is not None:
                                self.load_count_var.set(f"预计命中 {hit_count:,} 个实体细单元")
                            self.status_var.set("统计完成")
                    self._set_task_buttons("normal")
                elif kind == "done":
                    self.status_var.set("生成完成")
                    self._set_task_buttons("normal")
                elif kind == "solve_done":
                    _, seconds, returncode, output_dir = item
                    self._set_task_buttons("normal")
                    self.status_var.set("求解完成" if returncode == 0 else "求解失败")
                    self._append_log(f"求解耗时: {seconds:.3f} 秒；结果目录: {output_dir}")
                    if returncode:
                        messagebox.showerror("求解失败", f"DLUTFEM 返回码: {returncode}，请查看日志。")
                elif kind == "error":
                    self.status_var.set("失败")
                    self._set_task_buttons("normal")
                    self._append_log(item[1])
                    if len(item) > 2 and item[2]:
                        self._append_log(item[2])
                    messagebox.showerror("运行失败", item[1])
        except queue.Empty:
            pass
        self.after(100, self._poll_queue)


def main() -> None:
    app = VoxelGui()
    app.mainloop()


if __name__ == "__main__":
    main()
