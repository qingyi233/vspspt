"""Grid geometry shared by the GUI and command-line converter."""

from __future__ import annotations

import numpy as np


def validate_box(origin, size):
    origin = np.asarray(origin, dtype=float)
    size = np.asarray(size, dtype=float)
    if origin.shape != (3,) or size.shape != (3,):
        raise ValueError("包围盒位置和尺寸都必须包含 X/Y/Z 三个数值")
    if not np.isfinite(origin).all() or not np.isfinite(size).all() or not np.isfinite(origin + size).all():
        raise ValueError("包围盒位置和尺寸必须是有限数值")
    if np.any(size <= 0):
        raise ValueError("包围盒 X/Y/Z 尺寸必须全部大于 0")
    if np.any(origin + size <= origin):
        raise ValueError("包围盒尺寸过小，无法在给定位置表示")
    return tuple(float(x) for x in origin), tuple(float(x) for x in size)


def mesh_bounds(mesh_vertices):
    mesh = np.asarray(mesh_vertices)
    if mesh.ndim != 3 or mesh.shape[1:] != (3, 3) or not mesh.size or not np.isfinite(mesh).all():
        raise ValueError("模型必须包含有效三角面和有限坐标")
    return mesh.min(axis=(0, 2)), mesh.max(axis=(0, 2))


def resolve_grid(mesh_vertices, target_cells, longest_resolution, box_origin=None, box_size=None, fine=None):
    mesh_min, mesh_max = mesh_bounds(mesh_vertices)
    if (box_origin is None) != (box_size is None):
        raise ValueError("自定义包围盒必须同时指定位置和尺寸")
    if box_origin is None:
        origin, size = validate_box(mesh_min, mesh_max - mesh_min)
    else:
        origin, size = validate_box(box_origin, box_size)
        # The user requires containment, never clipping. Do not silently add tolerance.
        if np.any(np.asarray(origin) > mesh_min) or np.any(np.asarray(origin) + size < mesh_max):
            raise ValueError(
                "自定义包围盒未完全包含模型，已停止生成。\n"
                f"模型最小角点: {tuple(float(x) for x in mesh_min)}，最大角点: {tuple(float(x) for x in mesh_max)}\n"
                f"包围盒最小角点: {origin}，最大角点: {tuple(float(x) for x in np.asarray(origin) + size)}"
            )
    spans = np.asarray(size)
    ratios = spans / spans.max()
    if longest_resolution is not None:
        base = longest_resolution
    else:
        base = (target_cells / float(np.prod(ratios))) ** (1.0 / 3.0)
    dims = tuple(int(x) for x in np.maximum(2, np.rint(base * ratios)))
    # Keep an explicitly requested VGM box size exact without later boundary padding.
    if box_origin is not None and fine is not None:
        dims = tuple(((n + f - 1) // f) * f for n, f in zip(dims, fine))
    spacing = tuple(float(spans[axis] / dims[axis]) for axis in range(3))
    return dims, tuple(float(x) for x in origin), spacing


def add_box_arguments(parser):
    parser.add_argument("--box-origin", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="World-space minimum corner of the voxelization box; requires --box-size.")
    parser.add_argument("--box-size", nargs=3, type=float, metavar=("LX", "LY", "LZ"),
                        help="Positive box sizes. The box must fully contain the model; clipping is rejected.")
