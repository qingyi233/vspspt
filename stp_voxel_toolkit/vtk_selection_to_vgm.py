"""Convert voxel VTK plus two ParaView selections to a custom VGM file.

Workflow:
1. Generate occupancy and surface VTK with stp_to_voxel_vtk.py --surface.
2. Open the surface VTK in ParaView, select constraint faces, Extract Selection,
   and save it as a legacy .vtk file.
3. Do the same for load faces.
4. Run this script with the original occupancy VTK and the two selection VTKs.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

import numpy as np


VTK_DTYPE_BYTES = {
    "unsigned_char": 1,
    "uchar": 1,
    "char": 1,
    "int": 4,
    "unsigned_int": 4,
    "float": 4,
    "double": 8,
}

VTK_NUMPY_BE = {
    "unsigned_char": ">u1",
    "uchar": ">u1",
    "char": ">i1",
    "int": ">i4",
    "unsigned_int": ">u4",
    "float": ">f4",
    "double": ">f8",
}


@dataclass
class LegacyVtk:
    path: Path
    binary: bool
    dataset: str | None
    dimensions: tuple[int, int, int] | None
    origin: tuple[float, float, float] | None
    spacing: tuple[float, float, float] | None
    point_count: int | None
    cell_count: int | None
    cell_data: dict[str, np.ndarray]


def _readline(f: BinaryIO) -> bytes:
    line = f.readline()
    if not line:
        return b""
    return line


def _line_text(line: bytes) -> str:
    return line.decode("ascii", errors="ignore").strip()


def _skip_binary_payload(f: BinaryIO, byte_count: int) -> None:
    if byte_count > 0:
        f.seek(byte_count, 1)
    pos = f.tell()
    next_byte = f.read(1)
    if next_byte not in {b"", b"\n", b"\r"}:
        f.seek(pos)
    elif next_byte == b"\r":
        maybe_lf = f.read(1)
        if maybe_lf != b"\n":
            f.seek(f.tell() - 1)


def _read_ascii_numbers(f: BinaryIO, count: int, cast) -> np.ndarray:
    values = []
    while len(values) < count:
        line = _readline(f)
        if not line:
            break
        values.extend(_line_text(line).split())
    if len(values) < count:
        raise ValueError(f"Expected {count} ASCII values, got {len(values)}")
    return np.array([cast(x) for x in values[:count]])


def _read_data_array(f: BinaryIO, count: int, vtk_type: str, binary: bool) -> np.ndarray:
    vtk_type = vtk_type.lower()
    if vtk_type not in VTK_DTYPE_BYTES:
        raise ValueError(f"Unsupported VTK data type: {vtk_type}")
    if binary:
        byte_count = count * VTK_DTYPE_BYTES[vtk_type]
        data = f.read(byte_count)
        if len(data) != byte_count:
            raise ValueError(f"Expected {byte_count} bytes, got {len(data)}")
        arr = np.frombuffer(data, dtype=np.dtype(VTK_NUMPY_BE[vtk_type])).copy()
        _skip_binary_payload(f, 0)
        return arr
    cast = float if vtk_type in {"float", "double"} else int
    return _read_ascii_numbers(f, count, cast)


def read_legacy_vtk(path: Path) -> LegacyVtk:
    with path.open("rb") as f:
        version = _line_text(_readline(f))
        if not version.startswith("# vtk DataFile"):
            raise ValueError(f"{path} is not a legacy VTK file")
        _title = _readline(f)
        mode = _line_text(_readline(f)).upper()
        binary = mode == "BINARY"
        if mode not in {"BINARY", "ASCII"}:
            raise ValueError(f"Unsupported VTK mode in {path}: {mode}")

        dataset = None
        dimensions = None
        origin = None
        spacing = None
        point_count = None
        cell_count = None
        current_data_kind = None
        current_data_count = 0
        cell_data: dict[str, np.ndarray] = {}

        while True:
            line = _readline(f)
            if not line:
                break
            text = _line_text(line)
            if not text:
                continue
            parts = text.split()
            key = parts[0].upper()

            if key == "DATASET":
                dataset = parts[1].upper()
            elif key == "DIMENSIONS":
                dimensions = tuple(int(x) for x in parts[1:4])
                point_count = dimensions[0] * dimensions[1] * dimensions[2]
            elif key == "ORIGIN":
                origin = tuple(float(x) for x in parts[1:4])
            elif key == "SPACING":
                spacing = tuple(float(x) for x in parts[1:4])
            elif key == "POINTS":
                point_count = int(parts[1])
                vtk_type = parts[2]
                if binary:
                    _skip_binary_payload(f, point_count * 3 * VTK_DTYPE_BYTES[vtk_type.lower()])
                else:
                    _read_data_array(f, point_count * 3, vtk_type, binary=False)
            elif key in {"CELLS", "POLYGONS", "LINES", "VERTICES", "TRIANGLE_STRIPS"}:
                cell_count = int(parts[1])
                payload_count = int(parts[2])
                if binary:
                    _skip_binary_payload(f, payload_count * 4)
                else:
                    _read_ascii_numbers(f, payload_count, int)
            elif key == "CELL_TYPES":
                n = int(parts[1])
                if binary:
                    _skip_binary_payload(f, n * 4)
                else:
                    _read_ascii_numbers(f, n, int)
            elif key == "POINT_DATA":
                current_data_kind = "POINT_DATA"
                current_data_count = int(parts[1])
            elif key == "CELL_DATA":
                current_data_kind = "CELL_DATA"
                current_data_count = int(parts[1])
                cell_count = current_data_count
            elif key == "SCALARS":
                if len(parts) < 3:
                    raise ValueError(f"Malformed SCALARS line in {path}: {text}")
                name = parts[1]
                vtk_type = parts[2]
                components = int(parts[3]) if len(parts) > 3 else 1
                lookup = _line_text(_readline(f)).upper()
                if not lookup.startswith("LOOKUP_TABLE"):
                    raise ValueError(f"Expected LOOKUP_TABLE after SCALARS {name} in {path}")
                data = _read_data_array(f, current_data_count * components, vtk_type, binary)
                if components > 1:
                    data = data.reshape(current_data_count, components)
                if current_data_kind == "CELL_DATA":
                    cell_data[name] = data
            elif key == "FIELD":
                # Basic FIELD support for ParaView arrays such as vtkOriginalCellIds.
                field_count = int(parts[2])
                for _ in range(field_count):
                    meta = _line_text(_readline(f)).split()
                    if len(meta) < 4:
                        raise ValueError(f"Malformed FIELD metadata in {path}")
                    name, components_text, tuples_text, vtk_type = meta[:4]
                    components = int(components_text)
                    tuples = int(tuples_text)
                    data = _read_data_array(f, components * tuples, vtk_type, binary)
                    if components > 1:
                        data = data.reshape(tuples, components)
                    if current_data_kind == "CELL_DATA":
                        cell_data[name] = data

    return LegacyVtk(path, binary, dataset, dimensions, origin, spacing, point_count, cell_count, cell_data)


def read_occupancy(path: Path) -> tuple[np.ndarray, tuple[float, float, float], tuple[float, float, float]]:
    vtk = read_legacy_vtk(path)
    if vtk.dimensions is None or vtk.origin is None or vtk.spacing is None:
        raise ValueError("Occupancy VTK must contain DIMENSIONS, ORIGIN, and SPACING")
    if "occupancy" not in vtk.cell_data:
        raise ValueError("Occupancy VTK must contain CELL_DATA scalar named 'occupancy'")
    nx, ny, nz = (vtk.dimensions[0] - 1, vtk.dimensions[1] - 1, vtk.dimensions[2] - 1)
    expected = nx * ny * nz
    occ = vtk.cell_data["occupancy"]
    if occ.size != expected:
        raise ValueError(f"occupancy size mismatch: expected {expected}, got {occ.size}")
    return occ.astype(np.uint8).reshape((nx, ny, nz), order="F"), vtk.origin, vtk.spacing


def selected_voxels(path: Path, dims: tuple[int, int, int]) -> set[tuple[int, int, int]]:
    vtk = read_legacy_vtk(path)
    data = vtk.cell_data
    nx, ny, nz = dims
    if {"voxel_i", "voxel_j", "voxel_k"}.issubset(data):
        i = data["voxel_i"].astype(np.int64).ravel()
        j = data["voxel_j"].astype(np.int64).ravel()
        k = data["voxel_k"].astype(np.int64).ravel()
        return {
            (int(a), int(b), int(c))
            for a, b, c in zip(i, j, k)
            if 0 <= a < nx and 0 <= b < ny and 0 <= c < nz
        }
    if "voxel_eid" in data:
        eids = data["voxel_eid"].astype(np.int64).ravel()
    elif "vtkOriginalCellIds" in data:
        raise ValueError(
            f"{path.name} only has vtkOriginalCellIds. Regenerate the surface VTK with the updated tool "
            "so ParaView selections preserve voxel_i/voxel_j/voxel_k/voxel_eid."
        )
    else:
        available = ", ".join(sorted(data)) or "(none)"
        raise ValueError(
            f"{path.name} does not contain voxel_i/voxel_j/voxel_k or voxel_eid cell data. "
            f"Available cell data: {available}"
        )
    result = set()
    for eid in eids:
        zero = int(eid) - 1
        if zero < 0:
            continue
        k = zero // (nx * ny)
        rem = zero - k * nx * ny
        j = rem // nx
        i = rem - j * nx
        if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
            result.add((i, j, k))
    return result


def node_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return 1 + i + (nx + 1) * j + (nx + 1) * (ny + 1) * k


def coarse_node_id(i: int, j: int, k: int, cnx: int, cny: int) -> int:
    return 1 + i + (cnx + 1) * j + (cnx + 1) * (cny + 1) * k


def constraint_coarse_nodes(
    constraint_voxels: set[tuple[int, int, int]],
    fine: tuple[int, int, int],
    coarse_dims: tuple[int, int, int],
    mode: str,
) -> set[int]:
    fx, fy, fz = fine
    cnx, cny, cnz = coarse_dims
    coarse_cells = {
        (min(i // fx, cnx - 1), min(j // fy, cny - 1), min(k // fz, cnz - 1))
        for i, j, k in constraint_voxels
    }
    nodes: set[int] = set()
    if mode == "cell-corners":
        for ci, cj, ck in coarse_cells:
            for di in (0, 1):
                for dj in (0, 1):
                    for dk in (0, 1):
                        nodes.add(coarse_node_id(ci + di, cj + dj, ck + dk, cnx, cny))
    elif mode == "cell-min-corner":
        for ci, cj, ck in coarse_cells:
            nodes.add(coarse_node_id(ci, cj, ck, cnx, cny))
    else:
        raise ValueError(f"Unsupported SPC coarse-node mode: {mode}")
    return nodes


def bdf_field(value: object) -> str:
    text = str(value).strip()
    if text == "":
        return " " * 8
    if len(text) <= 8:
        return f"{text:>8}"
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(f"BDF small-field value is longer than 8 characters: {text!r}") from exc

    candidates = []
    for precision in range(6, 0, -1):
        candidates.append(f"{number:.{precision}E}".replace("E+0", "E+").replace("E-0", "E-"))
    candidates.extend(
        [
            f"{number:.7g}",
            f"{number:.6g}",
            f"{number:.5g}",
            f"{number:.4g}",
            f"{number:.3g}",
            f"{number:.2g}",
            f"{number:.1g}",
        ]
    )
    for candidate in candidates:
        candidate = candidate.replace("e+", "E+").replace("e-", "E-")
        if len(candidate) <= 8:
            return f"{candidate:>8}"
    raise ValueError(f"Cannot fit numeric value into 8-character BDF field: {text!r}")


def bdf_real(value: float) -> str:
    text = f"{value:.6g}"
    if "e" not in text.lower() and "." not in text:
        text += "."
    if len(text) <= 8:
        return text
    for precision in range(5, 0, -1):
        text = f"{value:.{precision}g}"
        if "e" not in text.lower() and "." not in text:
            text += "."
        if len(text) <= 8:
            return text
    text = f"{value:.1E}".replace("E+0", "E+").replace("E-0", "E-")
    if len(text) <= 8:
        return text
    raise ValueError(f"Cannot fit real value into 8-character BDF field: {value!r}")


def scalar_line(card: str, values: list[object], per_line: int = 8) -> list[str]:
    lines = []
    first = True
    for start in range(0, len(values), per_line):
        chunk = values[start : start + per_line]
        prefix = f"{card:<8}" if first else " " * 8
        fields = "".join(bdf_field(v) for v in chunk)
        lines.append((prefix + fields).rstrip())
        first = False
    return lines or [f"{card:<8}"]


def eid_vector_lines(card: str, eid: int, values: list[object], per_line: int = 7) -> list[str]:
    lines = []
    first_values = values[:per_line]
    lines.append((f"{card:<8}" + "".join(bdf_field(v) for v in [eid, *first_values])).rstrip())
    index = per_line
    while index < len(values):
        chunk = values[index : index + per_line]
        lines.append((" " * 16 + "".join(bdf_field(v) for v in chunk)).rstrip())
        index += per_line
    return lines


def material_type(block: np.ndarray) -> int:
    occupied = block.astype(bool, copy=False)
    count = int(occupied.sum())
    if count == 0:
        return 0
    if count == occupied.size:
        return 1
    return 2


def build_vgm(
    occupancy: np.ndarray,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    constraint_voxels: set[tuple[int, int, int]],
    load_voxels: set[tuple[int, int, int]],
    fine: tuple[int, int, int],
    spc_dof: str,
    force_sid: int,
    force_magnitude: float,
    force_direction: tuple[float, float, float],
    material_id: int,
    title: str,
    spc_node_mode: str = "cell-corners",
) -> list[str]:
    nx, ny, nz = occupancy.shape
    fx, fy, fz = fine
    cnx = math.ceil(nx / fx)
    cny = math.ceil(ny / fy)
    cnz = math.ceil(nz / fz)
    padded = np.zeros((cnx * fx, cny * fy, cnz * fz), dtype=np.uint8)
    padded[:nx, :ny, :nz] = occupancy
    fine_count = fx * fy * fz
    box_size = (nx * spacing[0], ny * spacing[1], nz * spacing[2])

    lines = [
        "$ VGM generated from voxel VTK and ParaView selections",
        f"$ ORIGIN {origin[0]:.17g} {origin[1]:.17g} {origin[2]:.17g}",
        "SOL 101",
        "CEND",
        f"TITLE = {title}",
        "SUBCASE 1",
        "   SUBTITLE=ParaView selections",
        "   SPC = 2",
        f"   LOAD = {force_sid + 20}",
        "   DISPLACEMENT(SORT1,REAL)=ALL",
        "   SPCFORCES(SORT1,REAL)=ALL",
        "   STRESS(SORT1,REAL,VONMISES,BILIN)=ALL",
    ]
    lines.extend(scalar_line("BOX", [bdf_real(box_size[0]), bdf_real(box_size[1]), bdf_real(box_size[2]), cnx, cny, cnz], 6))
    lines.extend(scalar_line("", [fx, fy, fz], 8))

    non_empty_eids: list[tuple[int, np.ndarray, int]] = []
    for ck in range(cnz):
        for cj in range(cny):
            for ci in range(cnx):
                eid = 1 + ci + cnx * cj + cnx * cny * ck
                block = padded[ci * fx : (ci + 1) * fx, cj * fy : (cj + 1) * fy, ck * fz : (ck + 1) * fz]
                mtype = material_type(block)
                lines.extend(scalar_line("COARSE", [eid, mtype, 0], 8))
                if mtype:
                    non_empty_eids.append((eid, block.copy(), mtype))

    lines.extend(scalar_line("MAT1", [material_id, "2.1E+11", "8.1E+10", "0.3", "7800."], 8))

    constraint_nodes = constraint_coarse_nodes(constraint_voxels, fine, (cnx, cny, cnz), spc_node_mode)
    lines.extend(scalar_line("SPCADD", [2, 1], 8))
    if constraint_nodes:
        lines.extend(scalar_line("SPC1", [1, spc_dof, *sorted(constraint_nodes)], 8))
    else:
        lines.append("$ No constraint selection was provided.")

    for eid, block, _mtype in non_empty_eids:
        values = block.ravel(order="F").astype(int).tolist()
        lines.extend(eid_vector_lines("MATB", eid, values))

    load_map: dict[int, np.ndarray] = {}
    for i, j, k in load_voxels:
        ci, cj, ck = i // fx, j // fy, k // fz
        eid = 1 + ci + cnx * cj + cnx * cny * ck
        local = (i % fx) + fx * (j % fy) + fx * fy * (k % fz)
        mask = load_map.setdefault(eid, np.zeros(fine_count, dtype=np.uint8))
        mask[local] = 1
    load_sid = force_sid + 20
    lines.extend(scalar_line("LOAD", [load_sid, "1.0", "1.0", force_sid], 8))
    header = [
        force_sid,
        0,
        bdf_real(force_magnitude),
        bdf_real(force_direction[0]),
        bdf_real(force_direction[1]),
        bdf_real(force_direction[2]),
        len(load_map),
    ]
    lines.extend(scalar_line("FORCE", header, 8))
    for eid in sorted(load_map):
        lines.extend(eid_vector_lines("", eid, load_map[eid].astype(int).tolist()))
    lines.append("ENDDATA")
    return lines


def write_vgm_stream(
    output_path: Path,
    occupancy: np.ndarray,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    constraint_voxels: set[tuple[int, int, int]],
    load_voxels: set[tuple[int, int, int]],
    fine: tuple[int, int, int],
    spc_dof: str,
    force_sid: int,
    force_magnitude: float,
    force_direction: tuple[float, float, float],
    material_id: int,
    title: str,
    progress=None,
    spc_node_mode: str = "cell-corners",
) -> None:
    nx, ny, nz = occupancy.shape
    fx, fy, fz = fine
    cnx = math.ceil(nx / fx)
    cny = math.ceil(ny / fy)
    cnz = math.ceil(nz / fz)
    padded = np.zeros((cnx * fx, cny * fy, cnz * fz), dtype=np.uint8)
    padded[:nx, :ny, :nz] = occupancy
    fine_count = fx * fy * fz
    box_size = (nx * spacing[0], ny * spacing[1], nz * spacing[2])
    total_coarse = cnx * cny * cnz

    constraint_nodes = constraint_coarse_nodes(constraint_voxels, fine, (cnx, cny, cnz), spc_node_mode)

    load_map: dict[int, np.ndarray] = {}
    for i, j, k in load_voxels:
        ci, cj, ck = i // fx, j // fy, k // fz
        eid = 1 + ci + cnx * cj + cnx * cny * ck
        local = (i % fx) + fx * (j % fy) + fx * fy * (k % fz)
        mask = load_map.setdefault(eid, np.zeros(fine_count, dtype=np.uint8))
        mask[local] = 1

    def write_lines(f, lines: list[str]) -> None:
        f.write(("\n".join(lines) + "\n").encode("ascii"))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as f:
        write_lines(
            f,
            [
                "$ VGM generated from voxel VTK and selected boundary conditions",
                f"$ ORIGIN {origin[0]:.17g} {origin[1]:.17g} {origin[2]:.17g}",
                "SOL 101",
                "CEND",
                f"TITLE = {title}",
                "SUBCASE 1",
                "   SUBTITLE=Generated boundary conditions",
                "   SPC = 2",
                f"   LOAD = {force_sid + 20}",
                "   DISPLACEMENT(SORT1,REAL)=ALL",
                "   SPCFORCES(SORT1,REAL)=ALL",
                "   STRESS(SORT1,REAL,VONMISES,BILIN)=ALL",
            ],
        )
        write_lines(f, scalar_line("BOX", [bdf_real(box_size[0]), bdf_real(box_size[1]), bdf_real(box_size[2]), cnx, cny, cnz], 6))
        write_lines(f, scalar_line("", [fx, fy, fz], 8))

        matb_records: list[tuple[int, np.ndarray]] = []
        processed = 0
        for ck in range(cnz):
            for cj in range(cny):
                for ci in range(cnx):
                    eid = 1 + ci + cnx * cj + cnx * cny * ck
                    block = padded[ci * fx : (ci + 1) * fx, cj * fy : (cj + 1) * fy, ck * fz : (ck + 1) * fz]
                    mtype = material_type(block)
                    write_lines(f, scalar_line("COARSE", [eid, mtype, 0], 8))
                    if mtype:
                        matb_records.append((eid, block.copy()))
                    processed += 1
                    if progress is not None and (processed % 100_000 == 0 or processed == total_coarse):
                        progress("COARSE", processed, total_coarse)

        write_lines(f, scalar_line("MAT1", [material_id, "2.1E+11", "8.1E+10", "0.3", "7800."], 8))
        write_lines(f, scalar_line("SPCADD", [2, 1], 8))
        if constraint_nodes:
            write_lines(f, scalar_line("SPC1", [1, spc_dof, *sorted(constraint_nodes)], 8))
        else:
            write_lines(f, ["$ No constraint selection was provided."])

        total_matb = len(matb_records)
        for index, (eid, block) in enumerate(matb_records, 1):
            write_lines(f, eid_vector_lines("MATB", eid, block.ravel(order="F").astype(int).tolist()))
            if progress is not None and (index % 10_000 == 0 or index == total_matb):
                progress("MATB", index, total_matb)

        load_sid = force_sid + 20
        write_lines(f, scalar_line("LOAD", [load_sid, "1.0", "1.0", force_sid], 8))
        header = [
            force_sid,
            0,
            bdf_real(force_magnitude),
            bdf_real(force_direction[0]),
            bdf_real(force_direction[1]),
            bdf_real(force_direction[2]),
            len(load_map),
        ]
        write_lines(f, scalar_line("FORCE", header, 8))
        total_load = len(load_map)
        for index, eid in enumerate(sorted(load_map), 1):
            write_lines(f, eid_vector_lines("", eid, load_map[eid].astype(int).tolist()))
            if progress is not None and (index % 5_000 == 0 or index == total_load):
                progress("FORCE", index, total_load)
        write_lines(f, ["ENDDATA"])


def vector3(text: str) -> tuple[float, float, float]:
    nums = [float(x) for x in re.split(r"[, ]+", text.strip()) if x]
    if len(nums) != 3:
        raise argparse.ArgumentTypeError("expected three numbers, e.g. 0,0,-1")
    return nums[0], nums[1], nums[2]


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert voxel VTK and two ParaView selection VTKs to VGM.")
    parser.add_argument("--occupancy", required=True, type=Path, help="Original *_occupancy_*.vtk file.")
    parser.add_argument("--constraint-selection", required=True, type=Path, help="ParaView Extract Selection VTK for constraints.")
    parser.add_argument("--load-selection", required=True, type=Path, help="ParaView Extract Selection VTK for loads.")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fine", nargs=3, type=int, default=(5, 5, 5), metavar=("NX", "NY", "NZ"))
    parser.add_argument("--spc-dof", default="123", help="SPC1 constrained DOF, e.g. 123.")
    parser.add_argument("--spc-node-mode", choices=("cell-corners", "cell-min-corner"), default="cell-corners")
    parser.add_argument("--force-sid", type=int, default=1)
    parser.add_argument("--force", type=float, default=1000.0, help="FORCE magnitude.")
    parser.add_argument("--direction", type=vector3, default=(0.0, 0.0, -1.0), help="Force direction, e.g. 0,0,-1.")
    parser.add_argument("--material-id", type=int, default=1)
    parser.add_argument("--title", default="Voxel VGM from ParaView selections")
    args = parser.parse_args()

    fine = tuple(int(x) for x in args.fine)
    if any(x <= 0 for x in fine):
        raise SystemExit("--fine values must be positive")

    occupancy, origin, spacing = read_occupancy(args.occupancy)
    dims = occupancy.shape
    constraint_voxels = selected_voxels(args.constraint_selection, dims)
    load_voxels = selected_voxels(args.load_selection, dims)
    def progress(stage, done, total):
        print(f"{stage}: {done:,}/{total:,}")

    write_vgm_stream(
        args.output,
        occupancy,
        origin,
        spacing,
        constraint_voxels,
        load_voxels,
        fine,
        args.spc_dof,
        args.force_sid,
        args.force,
        args.direction,
        args.material_id,
        args.title,
        progress=progress,
        spc_node_mode=args.spc_node_mode,
    )
    print(f"Wrote {args.output}")
    print(f"Grid: {dims[0]}x{dims[1]}x{dims[2]}, fine={fine[0]}x{fine[1]}x{fine[2]}")
    print(f"Constraint selected voxels: {len(constraint_voxels):,}")
    print(f"Load selected voxels: {len(load_voxels):,}")


if __name__ == "__main__":
    main()
