"""Mesh counts without constructing element connectivity or all node coordinates."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


def active_node_count(occupancy: np.ndarray) -> int:
    """Count the union of occupied hexahedra's corners using O(one plane) memory."""
    # Slice perpendicular to the longest axis to minimize temporary memory.
    cells = np.moveaxis(occupancy, int(np.argmax(occupancy.shape)), 0)
    a, b = cells.shape[1:]
    node_plane = np.zeros((a + 1, b + 1), dtype=bool)
    adjacent = np.zeros((a, b), dtype=bool)
    total = 0
    for index in range(cells.shape[0] + 1):
        adjacent.fill(False)
        if index:
            np.logical_or(adjacent, cells[index - 1], out=adjacent)
        if index < cells.shape[0]:
            np.logical_or(adjacent, cells[index], out=adjacent)
        node_plane.fill(False)
        node_plane[:-1, :-1] |= adjacent
        node_plane[1:, :-1] |= adjacent
        node_plane[:-1, 1:] |= adjacent
        node_plane[1:, 1:] |= adjacent
        total += int(np.count_nonzero(node_plane))
    return total


def mesh_statistics(occupancy, fine, origin, spacing, constraint_nodes=(), spc_dof="123") -> dict:
    if occupancy.ndim != 3 or any(n <= 0 for n in occupancy.shape):
        raise ValueError("统计需要非空三维占据场")
    coarse = occupancy
    for axis, step in enumerate(fine):
        coarse = np.maximum.reduceat(coarse, np.arange(0, coarse.shape[axis], step), axis=axis)
    coarse_shape = tuple(int(n) for n in coarse.shape)
    padded_shape = tuple(n * f for n, f in zip(coarse_shape, fine))
    fine_nodes = active_node_count(occupancy)
    coarse_nodes = active_node_count(coarse)
    return {
        "origin": [float(x) for x in origin],
        "spacing": [float(x) for x in spacing],
        "box_size": [n * float(s) for n, s in zip(padded_shape, spacing)],
        "input_dimensions": list(occupancy.shape),
        "fine_subdivisions": list(fine),
        "fine_grid_dimensions": list(padded_shape),
        "fine_grid_elements": math.prod(padded_shape),
        "fine_grid_nodes": math.prod(n + 1 for n in padded_shape),
        "fine_grid_dofs": math.prod(n + 1 for n in padded_shape) * 3,
        "solid_fine_elements": int(np.count_nonzero(occupancy)),
        "solid_fine_nodes": fine_nodes,
        "solid_fine_dofs": fine_nodes * 3,
        "coarse_grid_dimensions": list(coarse_shape),
        "coarse_grid_elements": math.prod(coarse_shape),
        "coarse_grid_nodes": math.prod(n + 1 for n in coarse_shape),
        "coarse_grid_dofs": math.prod(n + 1 for n in coarse_shape) * 3,
        "active_coarse_elements": int(np.count_nonzero(coarse)),
        "active_coarse_nodes": coarse_nodes,
        "active_coarse_dofs": coarse_nodes * 3,
        "constrained_coarse_nodes": len(constraint_nodes),
        "constrained_coarse_dofs": len(constraint_nodes) * len(set(spc_dof) & set("123")),
        "dof_definition": "3 translations per active node, before constraints or solver condensation",
    }


def save_statistics(output_path: Path, stats: dict) -> Path:
    path = output_path.with_suffix(".stats.json")
    path.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def format_statistics(stats: dict) -> str:
    return (
        f"实体细网格: {stats['solid_fine_elements']:,} 单元 / "
        f"{stats['solid_fine_nodes']:,} 节点 / {stats['solid_fine_dofs']:,} 自由度\n"
        f"有效粗网格: {stats['active_coarse_elements']:,} 单元 / "
        f"{stats['active_coarse_nodes']:,} 节点 / {stats['active_coarse_dofs']:,} 自由度\n"
        f"包围盒粗网格（含空单元）: {stats['coarse_grid_elements']:,} 单元 / "
        f"{stats['coarse_grid_nodes']:,} 节点\n"
        f"约束粗节点: {stats['constrained_coarse_nodes']:,}，约束位移自由度: "
        f"{stats['constrained_coarse_dofs']:,}\n"
        "自由度按每节点 3 个平移分量统计，未扣除约束或求解器凝聚。"
    )
