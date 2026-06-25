"""GUI for converting STEP/STP CAD files into voxelized VTK outputs.

Pipeline:
1. STEP/STP -> triangulated STL with Gmsh.
2. STL -> structured occupancy grid with pyvoxel/cvoxelise.
3. Occupancy grid -> legacy binary VTK STRUCTURED_POINTS.
4. Optional occupancy grid -> visible exterior surface VTK UNSTRUCTURED_GRID.
"""

from __future__ import annotations

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
from vtk_selection_to_vgm import read_occupancy, selected_voxels, write_vgm_stream
from run_dlutfem import DEFAULT_SOLVER_BIN, DEFAULT_WORK_ROOT, run_dlutfem, safe_job_name


ROOT = Path(__file__).resolve().parent
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
    out_dir: Path
    target_cells: int
    longest_resolution: int | None
    mesh_size: float | None
    ray: str
    reuse_stl: bool
    use_cpp: bool
    write_occupancy: bool
    write_surface: bool


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


def run_pipeline(config: RunConfig, q: queue.Queue) -> None:
    total_start = time.time()
    config.out_dir.mkdir(parents=True, exist_ok=True)

    input_path = config.step_path
    stl_path = config.out_dir / f"{input_path.stem}.stl"
    suffix = f"r{config.longest_resolution}" if config.longest_resolution else f"target{config.target_cells}"

    log_put(q, f"输入文件: {input_path}")
    log_put(q, f"输出目录: {config.out_dir}")

    start = time.time()
    if input_path.suffix.lower() == ".stl":
        stl_path = input_path
        stage_put(q, "STEP 转 STL", 0.0, "输入已是 STL，跳过")
    elif input_path.suffix.lower() in {".stp", ".step"} and config.reuse_stl and stl_path.exists() and stl_path.stat().st_size > 0:
        log_put(q, f"复用 STL: {stl_path}")
        stage_put(q, "STEP 转 STL", 0.0, "已复用")
    elif input_path.suffix.lower() in {".stp", ".step"}:
        triangles = step_to_stl(input_path, stl_path, config.mesh_size, q)
        stage_put(q, "STEP 转 STL", time.time() - start, f"{triangles:,} 个三角面, {stl_path.stat().st_size / 1048576:.2f} MB")
    else:
        raise ValueError("输入文件必须是 .stp、.step 或 .stl")

    start = time.time()
    mesh_vertices = pvoxel.Read_stl(str(stl_path))[0]
    if mesh_vertices is None or mesh_vertices.size == 0:
        raise RuntimeError(f"没有从 {stl_path} 读取到三角面")
    mins = tuple(float(mesh_vertices[:, axis, :].min()) for axis in range(3))
    maxs = tuple(float(mesh_vertices[:, axis, :].max()) for axis in range(3))
    spans = tuple(maxs[i] - mins[i] for i in range(3))
    nx, ny, nz = grid_for_mesh(mesh_vertices, config.target_cells, config.longest_resolution)
    total_cells = nx * ny * nz
    detail = f"{mesh_vertices.shape[0]:,} 个三角面, 网格 {nx}x{ny}x{nz}={total_cells:,}"
    stage_put(q, "读取 STL / 生成网格", time.time() - start, detail)
    log_put(q, f"包围盒最小值: {mins}")
    log_put(q, f"包围盒最大值: {maxs}")
    log_put(q, f"预计占据场内存: {total_cells / 1e6:.1f} MB")

    start = time.time()
    if config.use_cpp:
        vox_obj = voxel.Voxel(nx, ny, nz)
        grid = vox_obj.VOXELISE(mesh_vertices, ray=config.ray, parallel=False)
    else:
        vox_obj = pvoxel.Voxel(nx, ny, nz)
        grid = vox_obj.VOXELISE(mesh_vertices, ray=config.ray, parallel=False)
    grid = np.asarray(grid, dtype=np.uint8)
    if grid.ndim == 4:
        grid = grid[..., 0]
    occupied = int(grid.sum())
    stage_put(q, "体素化", time.time() - start, f"占据体素 {occupied:,}, 填充率 {occupied / total_cells:.8f}")

    spacing = tuple(spans[i] / max(1, [nx, ny, nz][i]) for i in range(3))
    size_tag = f"{nx}x{ny}x{nz}"
    occupancy_path = config.out_dir / f"{input_path.stem}_occupancy_{suffix}_{size_tag}.vtk"
    surface_path = config.out_dir / f"{input_path.stem}_surface_{suffix}_{size_tag}.vtk"

    if config.write_occupancy:
        start = time.time()
        write_occupancy_vtk(occupancy_path, grid, mins, spacing)
        stage_put(q, "写出占据场 VTK", time.time() - start, f"{occupancy_path.name}, {occupancy_path.stat().st_size / 1048576:.2f} MB")
        log_put(q, f"占据场 VTK: {occupancy_path}")
        q.put(("occupancy_path", str(occupancy_path)))

    if config.write_surface:
        start = time.time()
        face_count, point_count = write_surface_vtk(surface_path, grid, mins, spacing, q)
        detail = f"{face_count:,} 个四边形面, {point_count:,} 个点, {surface_path.stat().st_size / 1048576:.2f} MB"
        stage_put(q, "提取表面 VTK", time.time() - start, detail)
        log_put(q, f"可视表面 VTK: {surface_path}")

    stage_put(q, "总计", time.time() - total_start, "完成")
    q.put(("done", True))


class VoxelGui(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("STP/STL 体素化 VTK 转换工具")
        self.geometry("980x620")
        self.queue: queue.Queue = queue.Queue()
        self.worker: threading.Thread | None = None
        self.stage_rows: list[tuple[str, float, str]] = []

        default_step = ROOT.parent / "test1.stp"
        self.step_var = tk.StringVar(value=str(default_step if default_step.exists() else ROOT / "test1.stp"))
        self.out_var = tk.StringVar(value=str(ROOT / "output" / "gui_voxel_output"))
        self.target_var = tk.StringVar(value="100000000")
        self.resolution_var = tk.StringVar(value="")
        self.mesh_size_var = tk.StringVar(value="")
        self.ray_var = tk.StringVar(value="xyz")
        self.reuse_stl_var = tk.BooleanVar(value=True)
        self.use_cpp_var = tk.BooleanVar(value=True)
        self.write_occ_var = tk.BooleanVar(value=True)
        self.write_surface_var = tk.BooleanVar(value=True)
        self.vgm_occ_var = tk.StringVar(value="")
        self.vgm_constraint_var = tk.StringVar(value="")
        self.vgm_load_var = tk.StringVar(value="")
        self.vgm_output_var = tk.StringVar(value=str(ROOT / "output" / "model.vgm"))
        self.vgm_fine_x_var = tk.StringVar(value="5")
        self.vgm_fine_y_var = tk.StringVar(value="5")
        self.vgm_fine_z_var = tk.StringVar(value="5")
        self.vgm_spc_dof_var = tk.StringVar(value="123")
        self.vgm_force_var = tk.StringVar(value="1000")
        self.vgm_direction_var = tk.StringVar(value="0,0,-1")
        self.vgm_workers_var = tk.StringVar(value=str(max(1, min(8, (os.cpu_count() or 2) // 2))))
        self.vgm_spc_node_mode_var = tk.StringVar(value="命中粗单元8角点")
        self.vgm_button: ttk.Button | None = None
        self.solve_vgm_var = tk.StringVar(value="")
        self.solve_out_var = tk.StringVar(value=str(ROOT / "output" / "dlutfem_result"))
        self.solve_job_var = tk.StringVar(value="")
        self.solve_solver_bin_var = tk.StringVar(value=str(DEFAULT_SOLVER_BIN))
        self.solve_work_root_var = tk.StringVar(value=str(DEFAULT_WORK_ROOT))
        self.solve_refresh_bin_var = tk.BooleanVar(value=False)
        self.solve_button: ttk.Button | None = None

        self._build_ui()
        self.after(100, self._poll_queue)
        if IMPORT_ERROR is not None:
            messagebox.showerror("导入错误", f"缺少依赖或模块导入失败:\n{IMPORT_ERROR}")

    def _build_ui(self) -> None:
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)

        switch_frame = ttk.Frame(main)
        switch_frame.pack(fill="x")
        content_frame = ttk.Frame(main)
        content_frame.pack(fill="x", pady=(4, 0))
        voxel_tab = ttk.Frame(content_frame)
        vgm_tab = ttk.Frame(content_frame)
        solve_tab = ttk.Frame(content_frame)
        self._tab_buttons: dict[str, ttk.Button] = {}

        def show_tab(name: str) -> None:
            voxel_tab.pack_forget()
            vgm_tab.pack_forget()
            solve_tab.pack_forget()
            if name == "voxel":
                voxel_tab.pack(fill="x")
            elif name == "vgm":
                vgm_tab.pack(fill="x")
            else:
                solve_tab.pack(fill="x")
            self._tab_buttons["voxel"].configure(state="disabled" if name == "voxel" else "normal")
            self._tab_buttons["vgm"].configure(state="disabled" if name == "vgm" else "normal")
            self._tab_buttons["solve"].configure(state="disabled" if name == "solve" else "normal")

        self._tab_buttons["voxel"] = ttk.Button(switch_frame, text="体素化生成 VTK", command=lambda: show_tab("voxel"))
        self._tab_buttons["vgm"] = ttk.Button(switch_frame, text="选区生成 VGM", command=lambda: show_tab("vgm"))
        self._tab_buttons["solve"] = ttk.Button(switch_frame, text="VGM 有限元求解", command=lambda: show_tab("solve"))
        self._tab_buttons["voxel"].pack(side="left")
        self._tab_buttons["vgm"].pack(side="left", padx=(6, 0))
        self._tab_buttons["solve"].pack(side="left", padx=(6, 0))

        file_frame = ttk.LabelFrame(voxel_tab, text="输入 / 输出", padding=10)
        file_frame.pack(fill="x")

        ttk.Label(file_frame, text="STEP/STP/STL 文件").grid(row=0, column=0, sticky="w")
        ttk.Entry(file_frame, textvariable=self.step_var).grid(row=0, column=1, sticky="ew", padx=6)
        ttk.Button(file_frame, text="浏览", command=self._browse_step).grid(row=0, column=2)

        ttk.Label(file_frame, text="输出目录").grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(file_frame, textvariable=self.out_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(8, 0))
        ttk.Button(file_frame, text="浏览", command=self._browse_out).grid(row=1, column=2, pady=(8, 0))
        file_frame.columnconfigure(1, weight=1)

        opt_frame = ttk.LabelFrame(voxel_tab, text="体素化设置", padding=10)
        opt_frame.pack(fill="x", pady=(10, 0))

        ttk.Label(opt_frame, text="目标总体素数").grid(row=0, column=0, sticky="w")
        ttk.Entry(opt_frame, textvariable=self.target_var, width=16).grid(row=0, column=1, sticky="w", padx=6)
        ttk.Label(opt_frame, text="最长轴分辨率").grid(row=0, column=2, sticky="w", padx=(20, 0))
        ttk.Entry(opt_frame, textvariable=self.resolution_var, width=12).grid(row=0, column=3, sticky="w", padx=6)

        ttk.Label(opt_frame, text="Gmsh 网格尺寸").grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(opt_frame, textvariable=self.mesh_size_var, width=16).grid(row=1, column=1, sticky="w", padx=6, pady=(8, 0))
        ttk.Label(opt_frame, text="射线方向").grid(row=1, column=2, sticky="w", padx=(20, 0), pady=(8, 0))
        ttk.Combobox(opt_frame, textvariable=self.ray_var, values=("xyz", "x", "y", "z", "xy", "xz", "yz"), width=10, state="readonly").grid(row=1, column=3, sticky="w", padx=6, pady=(8, 0))

        check_frame = ttk.Frame(opt_frame)
        check_frame.grid(row=2, column=0, columnspan=4, sticky="w", pady=(10, 0))
        ttk.Checkbutton(check_frame, text="复用已有 STL", variable=self.reuse_stl_var).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(check_frame, text="使用 C++ 体素化", variable=self.use_cpp_var).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(check_frame, text="写出占据场 VTK", variable=self.write_occ_var).pack(side="left", padx=(0, 18))
        ttk.Checkbutton(check_frame, text="写出可视表面 VTK", variable=self.write_surface_var).pack(side="left")

        run_frame = ttk.Frame(voxel_tab)
        run_frame.pack(fill="x", pady=(10, 0))
        self.run_button = ttk.Button(run_frame, text="开始转换", command=self._start)
        self.run_button.pack(side="left")

        table_frame = ttk.LabelFrame(voxel_tab, text="耗时分析", padding=10)
        table_frame.pack(fill="x", pady=(10, 0))
        self.tree = ttk.Treeview(table_frame, columns=("stage", "seconds", "percent", "cumulative", "detail"), show="headings", height=6)
        self.tree.heading("stage", text="模块")
        self.tree.heading("seconds", text="耗时/秒")
        self.tree.heading("percent", text="占比")
        self.tree.heading("cumulative", text="累计/秒")
        self.tree.heading("detail", text="详情")
        self.tree.column("stage", width=150, anchor="w")
        self.tree.column("seconds", width=90, anchor="e")
        self.tree.column("percent", width=80, anchor="e")
        self.tree.column("cumulative", width=100, anchor="e")
        self.tree.column("detail", width=520, anchor="w")
        self.tree.pack(fill="x")

        vgm_frame = ttk.LabelFrame(vgm_tab, text="ParaView 选区生成 VGM", padding=10)
        vgm_frame.pack(fill="x")

        ttk.Label(vgm_frame, text="占据场 VTK").grid(row=0, column=0, sticky="w")
        ttk.Entry(vgm_frame, textvariable=self.vgm_occ_var).grid(row=0, column=1, sticky="ew", padx=6)
        ttk.Button(vgm_frame, text="浏览", command=self._browse_vgm_occ).grid(row=0, column=2)

        ttk.Label(vgm_frame, text="约束 VTK/STL").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(vgm_frame, textvariable=self.vgm_constraint_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(6, 0))
        ttk.Button(vgm_frame, text="浏览", command=self._browse_vgm_constraint).grid(row=1, column=2, pady=(6, 0))

        ttk.Label(vgm_frame, text="载荷 VTK/STL").grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(vgm_frame, textvariable=self.vgm_load_var).grid(row=2, column=1, sticky="ew", padx=6, pady=(6, 0))
        ttk.Button(vgm_frame, text="浏览", command=self._browse_vgm_load).grid(row=2, column=2, pady=(6, 0))

        ttk.Label(vgm_frame, text="输出 VGM").grid(row=3, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(vgm_frame, textvariable=self.vgm_output_var).grid(row=3, column=1, sticky="ew", padx=6, pady=(6, 0))
        ttk.Button(vgm_frame, text="浏览", command=self._browse_vgm_output).grid(row=3, column=2, pady=(6, 0))

        param_frame = ttk.Frame(vgm_frame)
        param_frame.grid(row=4, column=0, columnspan=3, sticky="w", pady=(8, 0))
        ttk.Label(param_frame, text="细分 X").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_fine_x_var, width=5).pack(side="left", padx=(6, 8))
        ttk.Label(param_frame, text="Y").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_fine_y_var, width=5).pack(side="left", padx=(6, 8))
        ttk.Label(param_frame, text="Z").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_fine_z_var, width=5).pack(side="left", padx=(6, 16))
        ttk.Label(param_frame, text="约束自由度").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_spc_dof_var, width=8).pack(side="left", padx=(6, 16))
        ttk.Label(param_frame, text="力大小").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_force_var, width=10).pack(side="left", padx=(6, 16))
        ttk.Label(param_frame, text="方向").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_direction_var, width=12).pack(side="left", padx=(6, 16))
        ttk.Label(param_frame, text="并行线程").pack(side="left")
        ttk.Entry(param_frame, textvariable=self.vgm_workers_var, width=5).pack(side="left", padx=(6, 16))
        self.vgm_button = ttk.Button(param_frame, text="生成 VGM", command=self._start_vgm)
        self.vgm_button.pack(side="left")
        vgm_frame.columnconfigure(1, weight=1)

        spc_mode_frame = ttk.Frame(vgm_frame)
        spc_mode_frame.grid(row=5, column=0, columnspan=3, sticky="w", pady=(8, 0))
        ttk.Label(spc_mode_frame, text="SPC粗节点映射").pack(side="left")
        ttk.Combobox(
            spc_mode_frame,
            textvariable=self.vgm_spc_node_mode_var,
            values=("命中粗单元8角点", "命中粗单元最小角点"),
            width=18,
            state="readonly",
        ).pack(side="left", padx=(6, 16))

        vgm_help = ttk.Label(vgm_frame, text="占据场必须选择 occupancy VTK；约束/载荷可以选择选区 VTK 或 UG 提取面的 STL。默认每个 VGM 粗单元包含 5x5x5 个细体素。", foreground="#555555")
        vgm_help.grid(row=6, column=0, columnspan=3, sticky="w", pady=(8, 0))

        solve_frame = ttk.LabelFrame(solve_tab, text="DLUTFEM 求解 VGM", padding=10)
        solve_frame.pack(fill="x")

        ttk.Label(solve_frame, text="输入 VGM").grid(row=0, column=0, sticky="w")
        ttk.Entry(solve_frame, textvariable=self.solve_vgm_var).grid(row=0, column=1, sticky="ew", padx=6)
        ttk.Button(solve_frame, text="浏览", command=self._browse_solve_vgm).grid(row=0, column=2)

        ttk.Label(solve_frame, text="结果输出目录").grid(row=1, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(solve_frame, textvariable=self.solve_out_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(6, 0))
        ttk.Button(solve_frame, text="浏览", command=self._browse_solve_out).grid(row=1, column=2, pady=(6, 0))

        ttk.Label(solve_frame, text="Job 名称").grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(solve_frame, textvariable=self.solve_job_var, width=24).grid(row=2, column=1, sticky="w", padx=6, pady=(6, 0))

        ttk.Label(solve_frame, text="求解器 bin 目录").grid(row=3, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(solve_frame, textvariable=self.solve_solver_bin_var).grid(row=3, column=1, sticky="ew", padx=6, pady=(6, 0))
        ttk.Button(solve_frame, text="浏览", command=self._browse_solve_bin).grid(row=3, column=2, pady=(6, 0))

        ttk.Label(solve_frame, text="英文临时目录").grid(row=4, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(solve_frame, textvariable=self.solve_work_root_var).grid(row=4, column=1, sticky="ew", padx=6, pady=(6, 0))
        ttk.Button(solve_frame, text="浏览", command=self._browse_solve_work_root).grid(row=4, column=2, pady=(6, 0))

        solve_options = ttk.Frame(solve_frame)
        solve_options.grid(row=5, column=0, columnspan=3, sticky="w", pady=(10, 0))
        ttk.Checkbutton(solve_options, text="重新复制求解器运行库", variable=self.solve_refresh_bin_var).pack(side="left", padx=(0, 16))
        self.solve_button = ttk.Button(solve_options, text="开始求解", command=self._start_solve)
        self.solve_button.pack(side="left")

        solve_help = ttk.Label(solve_frame, text="求解时会把 VGM 和 DLUTFEM 临时复制到纯英文路径运行，避免中文路径导致求解器崩溃。结果会回拷到上面的输出目录。", foreground="#555555")
        solve_help.grid(row=6, column=0, columnspan=3, sticky="w", pady=(8, 0))
        solve_frame.columnconfigure(1, weight=1)

        self.status_var = tk.StringVar(value="就绪")
        status_frame = ttk.Frame(main)
        status_frame.pack(fill="x", pady=(10, 0))
        ttk.Label(status_frame, text="状态").pack(side="left")
        ttk.Label(status_frame, textvariable=self.status_var).pack(side="left", padx=8)

        log_frame = ttk.LabelFrame(main, text="日志", padding=10)
        log_frame.pack(fill="both", expand=True, pady=(10, 0))
        self.log_text = tk.Text(log_frame, height=10, wrap="word")
        self.log_text.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=scrollbar.set)
        show_tab("voxel")

    def _browse_step(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("CAD/网格文件", "*.stp *.step *.stl"), ("STEP 文件", "*.stp *.step"), ("STL 文件", "*.stl"), ("所有文件", "*.*")])
        if path:
            self.step_var.set(path)

    def _browse_out(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.out_var.set(path)

    def _browse_vgm_occ(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("VTK 文件", "*.vtk"), ("所有文件", "*.*")])
        if path:
            self.vgm_occ_var.set(path)

    def _browse_vgm_constraint(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("VTK/STL 文件", "*.vtk *.stl"), ("VTK 文件", "*.vtk"), ("STL 文件", "*.stl"), ("所有文件", "*.*")])
        if path:
            self.vgm_constraint_var.set(path)

    def _browse_vgm_load(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("VTK/STL 文件", "*.vtk *.stl"), ("VTK 文件", "*.vtk"), ("STL 文件", "*.stl"), ("所有文件", "*.*")])
        if path:
            self.vgm_load_var.set(path)

    def _browse_vgm_output(self) -> None:
        path = filedialog.asksaveasfilename(defaultextension=".vgm", filetypes=[("VGM 文件", "*.vgm"), ("所有文件", "*.*")])
        if path:
            self.vgm_output_var.set(path)
            if not self.solve_vgm_var.get().strip():
                self.solve_vgm_var.set(path)
            if not self.solve_job_var.get().strip():
                self.solve_job_var.set(safe_job_name(Path(path).stem))

    def _browse_solve_vgm(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("VGM 文件", "*.vgm"), ("所有文件", "*.*")])
        if path:
            self.solve_vgm_var.set(path)
            if not self.solve_job_var.get().strip():
                self.solve_job_var.set(safe_job_name(Path(path).stem))
            if self.solve_out_var.get().strip() == str(ROOT / "output" / "dlutfem_result"):
                self.solve_out_var.set(str(Path(path).parent / "dlutfem_result"))

    def _browse_solve_out(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.solve_out_var.set(path)

    def _browse_solve_bin(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.solve_solver_bin_var.set(path)

    def _browse_solve_work_root(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.solve_work_root_var.set(path)

    def _parse_config(self) -> RunConfig:
        step = Path(self.step_var.get())
        if not step.exists():
            raise ValueError("输入文件不存在")
        if step.suffix.lower() not in {".stp", ".step", ".stl"}:
            raise ValueError("输入文件必须是 .stp、.step 或 .stl")
        target = parse_positive_int(self.target_var.get(), "目标总体素数")

        resolution_text = self.resolution_var.get().strip()
        resolution = parse_positive_int(resolution_text, "最长轴分辨率") if resolution_text else None
        if resolution is not None and resolution <= 1:
            raise ValueError("最长轴分辨率必须大于 1")

        mesh_text = self.mesh_size_var.get().strip()
        mesh_size = float(mesh_text) if mesh_text else None
        if mesh_size is not None and mesh_size <= 0:
            raise ValueError("Gmsh 网格尺寸必须大于 0")

        if not self.write_occ_var.get() and not self.write_surface_var.get():
            raise ValueError("至少选择一种输出类型")

        return RunConfig(
            step_path=step,
            out_dir=Path(self.out_var.get()),
            target_cells=target,
            longest_resolution=resolution,
            mesh_size=mesh_size,
            ray=self.ray_var.get(),
            reuse_stl=self.reuse_stl_var.get(),
            use_cpp=self.use_cpp_var.get(),
            write_occupancy=self.write_occ_var.get(),
            write_surface=self.write_surface_var.get(),
        )

    def _start(self) -> None:
        if self.worker and self.worker.is_alive():
            return
        if IMPORT_ERROR is not None:
            messagebox.showerror("导入错误", str(IMPORT_ERROR))
            return
        try:
            config = self._parse_config()
        except Exception as exc:
            messagebox.showerror("参数错误", str(exc))
            return

        self.tree.delete(*self.tree.get_children())
        self.stage_rows.clear()
        self.log_text.delete("1.0", "end")
        self._set_task_buttons("disabled")
        self.status_var.set("运行中...")
        self.worker = threading.Thread(target=self._run_worker, args=(config,), daemon=True)
        self.worker.start()

    def _parse_vgm_config(self) -> dict[str, object]:
        occupancy = Path(self.vgm_occ_var.get())
        constraint = Path(self.vgm_constraint_var.get())
        load = Path(self.vgm_load_var.get())
        output = Path(self.vgm_output_var.get())
        if not occupancy.exists():
            raise ValueError("占据场 VTK 不存在")
        if not constraint.exists():
            raise ValueError("约束 VTK/STL 不存在")
        if not load.exists():
            raise ValueError("载荷 VTK/STL 不存在")
        if occupancy.suffix.lower() != ".vtk":
            raise ValueError("占据场必须选择 *_occupancy_*.vtk")
        if "_surface_" in occupancy.name.lower():
            raise ValueError("占据场不能选择 surface VTK，请选择同一批生成的 *_occupancy_*.vtk")
        if constraint.suffix.lower() not in {".vtk", ".stl"}:
            raise ValueError("约束文件必须是 ParaView 选区 .vtk 或曲面 .stl")
        if load.suffix.lower() not in {".vtk", ".stl"}:
            raise ValueError("载荷文件必须是 ParaView 选区 .vtk 或曲面 .stl")
        if not output.name:
            raise ValueError("请设置输出 VGM 文件")

        fine = (
            parse_positive_int(self.vgm_fine_x_var.get(), "细分 X"),
            parse_positive_int(self.vgm_fine_y_var.get(), "细分 Y"),
            parse_positive_int(self.vgm_fine_z_var.get(), "细分 Z"),
        )
        if fine == (1, 1, 1):
            answer = messagebox.askyesno(
                "确认细分设置",
                "当前细分为 1 1 1，会把每个 VTK 体素都写成 VGM 粗体素。\n"
                "一亿体素会生成非常大的 VGM，写出时间也会很长。\n\n"
                "是否继续？\n\n建议改为 5 5 5。",
            )
            if not answer:
                raise ValueError("已取消生成。请把细分改为 5 5 5 后再试。")
        force = float(self.vgm_force_var.get().strip())
        direction = parse_vector3(self.vgm_direction_var.get(), "力方向")
        workers = parse_positive_int(self.vgm_workers_var.get(), "并行线程")
        spc_node_mode_label = self.vgm_spc_node_mode_var.get()
        spc_node_mode = {
            "命中粗单元8角点": "cell-corners",
            "命中粗单元最小角点": "cell-min-corner",
        }.get(spc_node_mode_label)
        if spc_node_mode is None:
            raise ValueError("SPC粗节点映射选项无效")
        spc_dof = self.vgm_spc_dof_var.get().strip()
        if not spc_dof:
            raise ValueError("约束自由度不能为空")
        return {
            "occupancy": occupancy,
            "constraint": constraint,
            "load": load,
            "output": output,
            "fine": fine,
            "force": force,
            "direction": direction,
            "workers": workers,
            "spc_node_mode": spc_node_mode,
            "spc_dof": spc_dof,
        }

    def _start_vgm(self) -> None:
        if self.worker and self.worker.is_alive():
            messagebox.showwarning("正在运行", "当前已有任务在运行，请等待完成。")
            return
        try:
            config = self._parse_vgm_config()
        except Exception as exc:
            messagebox.showerror("参数错误", str(exc))
            return

        self._set_task_buttons("disabled")
        self.status_var.set("正在生成 VGM...")
        self.worker = threading.Thread(target=self._run_vgm_worker, args=(config,), daemon=True)
        self.worker.start()

    def _parse_solve_config(self) -> dict[str, object]:
        vgm_path = Path(self.solve_vgm_var.get())
        output_dir = Path(self.solve_out_var.get())
        solver_bin = Path(self.solve_solver_bin_var.get())
        work_root = Path(self.solve_work_root_var.get())
        job = self.solve_job_var.get().strip() or safe_job_name(vgm_path.stem)
        if not vgm_path.exists():
            raise ValueError("输入 VGM 不存在")
        if vgm_path.suffix.lower() != ".vgm":
            raise ValueError("输入文件必须是 .vgm")
        if not (solver_bin / "DLUTFEM.exe").exists():
            raise ValueError("求解器 bin 目录中找不到 DLUTFEM.exe")
        try:
            str(work_root).encode("ascii")
        except UnicodeEncodeError as exc:
            raise ValueError("英文临时目录不能包含中文或其他非 ASCII 字符") from exc
        return {
            "vgm": vgm_path,
            "output_dir": output_dir,
            "solver_bin": solver_bin,
            "work_root": work_root,
            "job": job,
            "refresh_bin": self.solve_refresh_bin_var.get(),
        }

    def _start_solve(self) -> None:
        if self.worker and self.worker.is_alive():
            messagebox.showwarning("正在运行", "当前已有任务在运行，请等待完成。")
            return
        try:
            config = self._parse_solve_config()
        except Exception as exc:
            messagebox.showerror("参数错误", str(exc))
            return

        self.log_text.delete("1.0", "end")
        self._set_task_buttons("disabled")
        self.status_var.set("正在求解 VGM...")
        self.worker = threading.Thread(target=self._run_solve_worker, args=(config,), daemon=True)
        self.worker.start()

    def _run_vgm_worker(self, config: dict[str, object]) -> None:
        try:
            start = time.time()
            occupancy_path = config["occupancy"]
            constraint_path = config["constraint"]
            load_path = config["load"]
            output_path = config["output"]
            fine = config["fine"]
            workers = int(config["workers"])
            spc_node_mode = str(config["spc_node_mode"])

            self.queue.put(("status", "正在读取占据场 VTK..."))
            log_put(self.queue, f"读取占据场 VTK: {occupancy_path}")
            occupancy, origin, spacing = read_occupancy(occupancy_path)
            dims = occupancy.shape
            log_put(self.queue, f"体素尺寸: {dims[0]} x {dims[1]} x {dims[2]}")

            log_put(self.queue, f"读取约束选区 VTK: {constraint_path}")
            if constraint_path.suffix.lower() == ".stl":
                tolerance = (spacing[0] ** 2 + spacing[1] ** 2 + spacing[2] ** 2) ** 0.5
                log_put(self.queue, f"约束文件是 STL，按曲面距离映射体素，容差: {tolerance:.9g}，并行线程: {workers}")
                self.queue.put(("status", "正在映射约束 STL..."))

                def constraint_progress(done, total, selected):
                    log_put(self.queue, f"约束 STL 进度: {done:,}/{total:,} 个三角面，当前命中 {selected:,} 个体素")

                constraint_voxels = surface_stl_to_voxels(
                    constraint_path,
                    occupancy,
                    origin,
                    spacing,
                    tolerance,
                    progress=constraint_progress,
                    workers=workers,
                )
            else:
                self.queue.put(("status", "正在读取约束选区 VTK..."))
                constraint_voxels = selected_voxels(constraint_path, dims)
            log_put(self.queue, f"约束选中体素: {len(constraint_voxels):,}")

            log_put(self.queue, f"读取载荷选区 VTK: {load_path}")
            if load_path.suffix.lower() == ".stl":
                tolerance = (spacing[0] ** 2 + spacing[1] ** 2 + spacing[2] ** 2) ** 0.5
                log_put(self.queue, f"载荷文件是 STL，按曲面距离映射体素，容差: {tolerance:.9g}，并行线程: {workers}")
                self.queue.put(("status", "正在映射载荷 STL..."))

                def load_progress(done, total, selected):
                    log_put(self.queue, f"载荷 STL 进度: {done:,}/{total:,} 个三角面，当前命中 {selected:,} 个体素")

                load_voxels = surface_stl_to_voxels(
                    load_path,
                    occupancy,
                    origin,
                    spacing,
                    tolerance,
                    progress=load_progress,
                    workers=workers,
                )
            else:
                self.queue.put(("status", "正在读取载荷选区 VTK..."))
                load_voxels = selected_voxels(load_path, dims)
            log_put(self.queue, f"载荷选中体素: {len(load_voxels):,}")

            self.queue.put(("status", "正在写出 VGM..."))
            log_put(self.queue, "开始流式写出 VGM...")

            def write_progress(stage, done, total):
                log_put(self.queue, f"写出 {stage}: {done:,}/{total:,}")

            write_vgm_stream(
                output_path,
                occupancy,
                origin,
                spacing,
                constraint_voxels,
                load_voxels,
                fine,
                str(config["spc_dof"]),
                1,
                float(config["force"]),
                config["direction"],
                1,
                "Voxel VGM from selected STL surfaces",
                progress=write_progress,
                spc_node_mode=spc_node_mode,
            )
            seconds = time.time() - start
            log_put(self.queue, f"VGM 输出: {output_path}")
            log_put(self.queue, f"VGM 文件大小: {output_path.stat().st_size / 1048576:.2f} MB")
            if not self.solve_vgm_var.get().strip():
                self.queue.put(("solve_vgm_path", str(output_path)))
            self.queue.put(("vgm_done", seconds, str(output_path)))
        except Exception:
            self.queue.put(("error", traceback.format_exc(), None))

    def _run_solve_worker(self, config: dict[str, object]) -> None:
        try:
            start = time.time()

            def log_callback(message: str) -> None:
                if message:
                    log_put(self.queue, message)

            returncode = run_dlutfem(
                config["vgm"],
                config["solver_bin"],
                config["work_root"],
                config["output_dir"],
                config["job"],
                bool(config["refresh_bin"]),
                log=log_callback,
            )
            seconds = time.time() - start
            self.queue.put(("solve_done", seconds, returncode, str(config["output_dir"])))
        except Exception:
            self.queue.put(("error", traceback.format_exc(), None))

    def _run_worker(self, config: RunConfig) -> None:
        try:
            run_pipeline(config, self.queue)
        except StepMeshingError as exc:
            self.queue.put(("error", str(exc), traceback.format_exc()))
        except Exception:
            self.queue.put(("error", traceback.format_exc(), None))

    def _set_task_buttons(self, state: str) -> None:
        self.run_button.configure(state=state)
        if self.vgm_button is not None:
            self.vgm_button.configure(state=state)
        if self.solve_button is not None:
            self.solve_button.configure(state=state)

    def _poll_queue(self) -> None:
        try:
            while True:
                item = self.queue.get_nowait()
                kind = item[0]
                if kind == "log":
                    self.log_text.insert("end", item[1] + "\n")
                    self.log_text.see("end")
                elif kind == "status":
                    self.status_var.set(item[1])
                elif kind == "stage":
                    _kind, stage, seconds, detail = item
                    self.stage_rows.append((stage, seconds, detail))
                    self._refresh_time_table(final=False)
                elif kind == "done":
                    self._refresh_time_table(final=True)
                    self.status_var.set("已完成")
                    self._set_task_buttons("normal")
                    messagebox.showinfo("完成", "转换完成。")
                elif kind == "vgm_done":
                    _kind, seconds, output_path = item
                    self.status_var.set("VGM 已完成")
                    self._set_task_buttons("normal")
                    self.log_text.insert("end", f"生成 VGM 耗时: {seconds:.3f} 秒\n")
                    self.log_text.see("end")
                    messagebox.showinfo("完成", f"VGM 生成完成:\n{output_path}")
                elif kind == "occupancy_path":
                    if not self.vgm_occ_var.get().strip():
                        self.vgm_occ_var.set(item[1])
                elif kind == "solve_vgm_path":
                    self.solve_vgm_var.set(item[1])
                    if not self.solve_job_var.get().strip():
                        self.solve_job_var.set(safe_job_name(Path(item[1]).stem))
                elif kind == "solve_done":
                    _kind, seconds, returncode, output_dir = item
                    self._set_task_buttons("normal")
                    if returncode == 0:
                        self.status_var.set("求解已完成")
                        self.log_text.insert("end", f"求解耗时: {seconds:.3f} 秒\n")
                        self.log_text.see("end")
                        messagebox.showinfo("完成", f"求解完成，结果目录:\n{output_dir}")
                    else:
                        self.status_var.set("求解失败")
                        self.log_text.insert("end", f"求解失败，返回码: {returncode}，耗时: {seconds:.3f} 秒\n")
                        self.log_text.see("end")
                        messagebox.showerror("求解失败", f"DLUTFEM 返回码: {returncode}\n请查看日志输出。")
                elif kind == "error":
                    self.status_var.set("失败")
                    self._set_task_buttons("normal")
                    self.log_text.insert("end", item[1] + "\n")
                    if len(item) > 2 and item[2]:
                        self.log_text.insert("end", "\n详细调试信息:\n" + item[2] + "\n")
                    self.log_text.see("end")
                    messagebox.showerror("转换失败", item[1])
        except queue.Empty:
            pass
        self.after(100, self._poll_queue)

    def _refresh_time_table(self, final: bool) -> None:
        self.tree.delete(*self.tree.get_children())
        non_total = [row for row in self.stage_rows if row[0] != "总计"]
        total_rows = [row for row in self.stage_rows if row[0] == "总计"]
        total_seconds = total_rows[-1][1] if total_rows else sum(row[1] for row in non_total)
        cumulative = 0.0

        for stage, seconds, detail in self.stage_rows:
            if stage != "总计":
                cumulative += seconds
            percent = (seconds / total_seconds * 100.0) if final and total_seconds > 0 else 0.0
            cumulative_text = f"{cumulative:.3f}" if stage != "总计" else f"{total_seconds:.3f}"
            percent_text = f"{percent:.1f}%" if final else "-"
            self.tree.insert(
                "",
                "end",
                values=(stage, f"{seconds:.3f}", percent_text, cumulative_text, detail),
            )


def main() -> None:
    app = VoxelGui()
    app.mainloop()


if __name__ == "__main__":
    main()
