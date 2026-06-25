"""GUI for converting STEP/STP CAD files into voxelized VTK outputs.

Pipeline:
1. STEP/STP -> triangulated STL with Gmsh.
2. STL -> structured occupancy grid with pyvoxel/cvoxelise.
3. Occupancy grid -> legacy binary VTK STRUCTURED_POINTS.
4. Optional occupancy grid -> visible exterior surface VTK UNSTRUCTURED_GRID.
"""

from __future__ import annotations

import queue
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
    log_put(q, "表面面片: " + ", ".join(f"{d}={c:,}" for d, c in counts))

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
        self.geometry("980x720")
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

        self._build_ui()
        self.after(100, self._poll_queue)
        if IMPORT_ERROR is not None:
            messagebox.showerror("导入错误", f"缺少依赖或模块导入失败:\n{IMPORT_ERROR}")

    def _build_ui(self) -> None:
        main = ttk.Frame(self, padding=12)
        main.pack(fill="both", expand=True)

        file_frame = ttk.LabelFrame(main, text="输入 / 输出", padding=10)
        file_frame.pack(fill="x")

        ttk.Label(file_frame, text="STEP/STP/STL 文件").grid(row=0, column=0, sticky="w")
        ttk.Entry(file_frame, textvariable=self.step_var).grid(row=0, column=1, sticky="ew", padx=6)
        ttk.Button(file_frame, text="浏览", command=self._browse_step).grid(row=0, column=2)

        ttk.Label(file_frame, text="输出目录").grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Entry(file_frame, textvariable=self.out_var).grid(row=1, column=1, sticky="ew", padx=6, pady=(8, 0))
        ttk.Button(file_frame, text="浏览", command=self._browse_out).grid(row=1, column=2, pady=(8, 0))
        file_frame.columnconfigure(1, weight=1)

        opt_frame = ttk.LabelFrame(main, text="体素化设置", padding=10)
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

        run_frame = ttk.Frame(main)
        run_frame.pack(fill="x", pady=(10, 0))
        self.run_button = ttk.Button(run_frame, text="开始转换", command=self._start)
        self.run_button.pack(side="left")
        self.status_var = tk.StringVar(value="就绪")
        ttk.Label(run_frame, textvariable=self.status_var).pack(side="left", padx=12)

        table_frame = ttk.LabelFrame(main, text="耗时分析", padding=10)
        table_frame.pack(fill="x", pady=(10, 0))
        self.tree = ttk.Treeview(table_frame, columns=("stage", "seconds", "percent", "cumulative", "detail"), show="headings", height=8)
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

        log_frame = ttk.LabelFrame(main, text="日志", padding=10)
        log_frame.pack(fill="both", expand=True, pady=(10, 0))
        self.log_text = tk.Text(log_frame, height=16, wrap="word")
        self.log_text.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=scrollbar.set)

    def _browse_step(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("CAD/网格文件", "*.stp *.step *.stl"), ("STEP 文件", "*.stp *.step"), ("STL 文件", "*.stl"), ("所有文件", "*.*")])
        if path:
            self.step_var.set(path)

    def _browse_out(self) -> None:
        path = filedialog.askdirectory()
        if path:
            self.out_var.set(path)

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
        self.run_button.configure(state="disabled")
        self.status_var.set("运行中...")
        self.worker = threading.Thread(target=self._run_worker, args=(config,), daemon=True)
        self.worker.start()

    def _run_worker(self, config: RunConfig) -> None:
        try:
            run_pipeline(config, self.queue)
        except StepMeshingError as exc:
            self.queue.put(("error", str(exc), traceback.format_exc()))
        except Exception:
            self.queue.put(("error", traceback.format_exc(), None))

    def _poll_queue(self) -> None:
        try:
            while True:
                item = self.queue.get_nowait()
                kind = item[0]
                if kind == "log":
                    self.log_text.insert("end", item[1] + "\n")
                    self.log_text.see("end")
                elif kind == "stage":
                    _kind, stage, seconds, detail = item
                    self.stage_rows.append((stage, seconds, detail))
                    self._refresh_time_table(final=False)
                elif kind == "done":
                    self._refresh_time_table(final=True)
                    self.status_var.set("已完成")
                    self.run_button.configure(state="normal")
                    messagebox.showinfo("完成", "转换完成。")
                elif kind == "error":
                    self.status_var.set("失败")
                    self.run_button.configure(state="normal")
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
