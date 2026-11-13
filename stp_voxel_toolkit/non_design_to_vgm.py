"""Select existing material cells in an analysis VGM using a closed STL region."""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np

from stl_surface_to_vgm import read_ascii_or_binary_stl
from vtk_selection_to_vgm import eid_vector_lines, scalar_line


@dataclass(frozen=True)
class AnalysisGrid:
    origin: tuple[float, float, float]
    lengths: tuple[float, float, float]
    coarse: tuple[int, int, int]
    fine: tuple[int, int, int]
    box_lines: tuple[str, str]


def _number(field):
    field = field.strip().upper().replace("D", "E")
    try:
        return float(field)
    except ValueError:
        # Fixed-field VGM also allows an exponent with no E, such as 1.25-3.
        return float(re.sub(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$", r"\1E\2", field))


def _fields(line, start=8):
    return [line[i:i + 8].strip() for i in range(start, min(len(line), 72), 8)]


def _cards(path):
    """Read fixed-field grid cards, including blank and + continuation fields."""
    name, lines = None, []
    with path.open("r", encoding="utf-8-sig", errors="replace") as source:
        for raw in source:
            line = raw.partition("$")[0].rstrip()
            if not line:
                comment = raw.lstrip()
                if comment.upper().startswith("$ ORIGIN "):
                    yield "ORIGIN", [comment]
                continue
            key = line[:8].strip().upper()
            if not key or key.startswith("+"):
                if name is not None:
                    lines.append(line)
                continue
            if name is not None:
                yield name, lines
            if key in {"BOX", "COARSE", "MATB"}:
                name, lines = key, [line]
            else:
                name, lines = None, []
        if name is not None:
            yield name, lines


def read_analysis_grid(path: Path, origin=None) -> AnalysisGrid:
    stored_origin, box = None, None
    for name, lines in _cards(path):
        if name == "ORIGIN":
            stored_origin = tuple(map(_number, lines[0].split()[2:]))
        elif name == "BOX":
            box = lines
        elif name in {"COARSE", "MATB"}:
            break
    if box is None or len(box) != 2:
        raise ValueError("Analysis VGM must contain a two-line fixed-field BOX card")
    fields = _fields(box[0])
    lengths = tuple(map(_number, fields[:3]))
    coarse = tuple(map(int, fields[3:6]))
    fine = tuple(map(int, _fields(box[1])[:3]))
    if (len(lengths) != 3 or len(coarse) != 3 or len(fine) != 3
            or any(not math.isfinite(value) or value <= 0 for value in lengths)
            or any(value <= 0 for value in coarse + fine)):
        raise ValueError("Invalid lengths or subdivisions in analysis VGM BOX card")
    if stored_origin is not None:
        if origin is not None and tuple(origin) != stored_origin:
            raise ValueError("Entered origin differs from $ ORIGIN in the analysis VGM")
        origin = stored_origin
    elif origin is not None:
        origin = tuple(origin)
    if origin is None or len(origin) != 3 or not all(map(math.isfinite, origin)):
        raise ValueError("Analysis VGM needs a $ ORIGIN comment or an explicit origin X Y Z")
    return AnalysisGrid(origin, lengths, coarse, fine, tuple(box))


def _region_mask(grid, geometry):
    runtime = str(Path(__file__).resolve().parent / "pyvoxel_runtime")
    if runtime not in sys.path:
        sys.path.insert(0, runtime)
    from voxel import Voxel

    triangles = read_ascii_or_binary_stl(geometry)
    if not triangles.size or not np.isfinite(triangles).all():
        raise ValueError("Non-design STL must contain finite triangles")
    _, indices = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    faces = indices.reshape(-1, 3)
    edges = np.sort(np.stack((faces, np.roll(faces, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
    _, edge_counts = np.unique(edges, axis=0, return_counts=True)
    if np.any(edges[:, 0] == edges[:, 1]) or np.any(edge_counts != 2):
        raise ValueError("Non-design STL must be a closed triangle mesh")

    mesh = triangles.transpose(0, 2, 1) - np.asarray(grid.origin)[None, :, None]
    low, high = mesh.min(axis=(0, 2)), mesh.max(axis=(0, 2))
    if np.any(high <= 0) or np.any(low >= grid.lengths):
        raise ValueError("Non-design STL does not overlap the analysis BOX")
    coarse_spacing = np.asarray(grid.lengths) / grid.coarse
    first = np.floor(low / coarse_spacing).astype(np.intp)
    end = np.ceil(high / coarse_spacing).astype(np.intp)
    first -= first * coarse_spacing > low
    end += end * coarse_spacing < high
    clipped_first = np.maximum(first, 0)
    clipped_end = np.minimum(end, grid.coarse)
    if np.any(clipped_first >= clipped_end):
        raise ValueError("Non-design STL does not overlap the analysis BOX")

    # Keep the crop aligned to the analysis cells, even when the STL extends past BOX.
    dims = tuple(map(int, (end - first) * grid.fine))
    bounds = np.stack((first * coarse_spacing, end * coarse_spacing))
    full = np.asarray(Voxel(*dims).VOXELISE(mesh, ray="xyz", parallel=False, bounds=bounds)).reshape(dims)
    candidates = int(np.count_nonzero(full))
    crop = tuple(slice(int((a - start) * f), int((b - start) * f))
                 for a, b, start, f in zip(clipped_first, clipped_end, first, grid.fine))
    return full[crop].astype(bool), clipped_first, candidates


def _intersect_analysis(path, grid, mask, first):
    # Each block view keeps the original X-fastest numbering inside its coarse cell.
    counts = tuple(n // f for n, f in zip(mask.shape, grid.fine))
    fx, fy, fz = grid.fine
    blocks = mask.reshape(counts[0], fx, counts[1], fy, counts[2], fz).transpose(0, 2, 4, 1, 3, 5)
    active = blocks.any(axis=(3, 4, 5))
    cell_types = np.full(counts, -1, dtype=np.int8)
    mixed_read = np.zeros(counts, dtype=bool)
    pending = {}
    coarse_count = math.prod(grid.coarse)
    nx, ny, _ = grid.coarse

    def apply_matb(index, density, eid):
        if len(density) != fx * fy * fz:
            raise ValueError(f"MATB {eid} has {len(density)} values, expected {fx * fy * fz}")
        if mixed_read[index]:
            raise ValueError(f"Duplicate MATB {eid}")
        blocks[index] &= np.array(density, dtype=bool).reshape(grid.fine, order="F")
        mixed_read[index] = True

    for name, lines in _cards(path):
        if name not in {"COARSE", "MATB"}:
            continue
        fields = _fields(lines[0])
        eid = int(fields[0])
        if not 1 <= eid <= coarse_count:
            raise ValueError(f"{name} element {eid} is outside analysis BOX numbering")
        z, plane = divmod(eid - 1, nx * ny)
        y, x = divmod(plane, nx)
        index = (x - int(first[0]), y - int(first[1]), z - int(first[2]))
        if any(i < 0 or i >= n for i, n in zip(index, counts)) or not active[index]:
            continue
        if name == "COARSE":
            cell_type = int(fields[1])
            if cell_types[index] != -1 or cell_type not in (0, 1, 2):
                raise ValueError(f"Duplicate or unsupported COARSE {eid}")
            cell_types[index] = cell_type
            if index in pending:
                density = pending.pop(index)
                if cell_type == 2:
                    apply_matb(index, density, eid)
        elif cell_types[index] in (-1, 2):
            entries = fields[1:] + [field for line in lines[1:] for field in _fields(line, 16)]
            density = [_number(field) > 0 if field else False for field in entries]
            if cell_types[index] == 2:
                apply_matb(index, density, eid)
            else:
                if index in pending:
                    raise ValueError(f"Duplicate MATB {eid}")
                pending[index] = density
    if pending or np.any((cell_types == 2) & ~mixed_read):
        raise ValueError("Selected mixed coarse cells require both COARSE and complete MATB cards")
    blocks[cell_types <= 0] = False
    return blocks


def _write_selection(output, grid, blocks, first):
    # Write a complete grid descriptor with a sparse material mask at original EIDs.
    nx, ny, _ = grid.coarse
    counts = blocks.sum(axis=(3, 4, 5))
    full_count = math.prod(grid.fine)
    types = {}
    for local in np.argwhere(counts):
        x, y, z = local + first
        eid = int(1 + x + nx * (y + ny * z))
        index = tuple(local)
        types[eid] = (1 if counts[index] == full_count else 2, index)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(mode="w", encoding="ascii", newline="\n", dir=output.parent,
                                suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            origin = " ".join(f"{x:.17g}" for x in grid.origin)
            stream.write(f"$ Non-design mask selected from analysis material\n$ ORIGIN {origin}\n"
                         "SOL 101\nCEND\nTITLE = Non-design domain\nSUBCASE 1\n")
            stream.write("\n".join(grid.box_lines) + "\n")
            for eid in range(1, math.prod(grid.coarse) + 1):
                cell_type = types[eid][0] if eid in types else 0
                stream.write(f"COARSE  {eid:8d}{cell_type:8d}{0:8d}\n")
            stream.write("\n".join(scalar_line("MAT1", [1, "2.1E+11", "8.1E+10", "0.3", "7800."])) + "\n")
            for eid in sorted(types):
                _, index = types[eid]
                values = blocks[index].ravel(order="F").astype(int).tolist()
                stream.write("\n".join(eid_vector_lines("MATB", eid, values)) + "\n")
            stream.write("ENDDATA\n")
        temporary.replace(output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return int(counts.sum()), len(types)


def select_non_design(analysis: Path, geometry: Path, output: Path, *, origin=None, progress=None):
    for input_path in (analysis, geometry):
        if output.resolve() == input_path.resolve() or (output.exists() and output.samefile(input_path)):
            raise ValueError("Non-design output must be separate from the analysis VGM and geometry STL")
    grid = read_analysis_grid(analysis, origin)
    if progress:
        progress("按分析网格计算非设计域几何")
    mask, first, candidates = _region_mask(grid, geometry)
    if progress:
        progress("与分析 VGM 实体细单元取交集")
    blocks = _intersect_analysis(analysis, grid, mask, first)
    if not blocks.any():
        raise ValueError("Non-design STL did not select any solid fine cells in the analysis VGM")
    solid, coarse = _write_selection(output, grid, blocks, first)
    return {"candidate_cells": candidates, "selected_cells": solid,
            "excluded_cells": candidates - solid, "selected_coarse_cells": coarse}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True, help="Analysis VGM defining the grid and material")
    parser.add_argument("--geometry", type=Path, required=True, help="Closed non-design region STL")
    parser.add_argument("--output", type=Path, required=True, help="Output non-design VGM mask")
    parser.add_argument("--origin", nargs=3, type=float, help="World-space origin when the VGM lacks $ ORIGIN")
    args = parser.parse_args()
    try:
        result = select_non_design(args.analysis, args.geometry, args.output, origin=args.origin)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"{exc}\n")
    print(f"Selected {result['selected_cells']:,} solid fine cells in {result['selected_coarse_cells']:,} coarse cells. "
          f"Excluded {result['excluded_cells']:,} candidates outside analysis material.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
