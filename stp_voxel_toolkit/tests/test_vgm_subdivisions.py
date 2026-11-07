"""Run with: python -m unittest discover -s stp_voxel_toolkit/tests -v"""

import argparse
import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

TOOLKIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLKIT))
import vtk_selection_to_vgm as vgm


def cards(text):
    result = []
    for line in text.splitlines():
        if line.startswith("$") or not line.strip():
            continue
        if line[:8].strip():
            result.append((line[:8].strip(), []))
        result[-1][1].append(line)
    return result


def fields(lines):
    return [line[i:i + 8].strip() for line in lines for i in range(8, len(line), 8)
            if line[i:i + 8].strip()]


class SubdivisionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.occ = np.ones((21, 11, 12), dtype=np.uint8)
        self.occ[10:20, :10, :10] = 0  # An empty coarse cell at fine=10.
        self.occ[2, 3, 4] = 0  # Mixed material and asymmetric local index.
        self.loads = {(9, 9, 9), (20, 10, 11), (0, 1, 2)}

    def check_geometry_and_loads(self, text, fine, occupancy=None, loads=None):
        occupancy = self.occ if occupancy is None else occupancy
        loads = self.loads if loads is None else loads
        parsed = cards(text)
        box = fields(next(lines for name, lines in parsed if name == "BOX"))
        np.testing.assert_array_equal(np.array(box[6:9], int), fine)
        coarse = np.array(box[3:6], int)
        shape = coarse * fine
        np.testing.assert_allclose(np.array(box[:3], float) / shape, (0.25, 0.5, 0.75))
        rebuilt = np.zeros(shape, dtype=np.uint8)
        decoded_loads = set()
        size = int(np.prod(fine))

        def location(eid, local):
            ci = (eid - 1) % coarse[0]
            cj = (eid - 1) // coarse[0] % coarse[1]
            ck = (eid - 1) // (coarse[0] * coarse[1])
            return (ci * fine[0] + local % fine[0],
                    cj * fine[1] + local // fine[0] % fine[1],
                    ck * fine[2] + local // (fine[0] * fine[1]))

        for name, lines in parsed:
            if name == "MATB":
                values = list(map(int, fields(lines)))
                self.assertEqual(len(values), size + 1)
                for local, value in enumerate(values[1:]):
                    rebuilt[location(values[0], local)] = value
            elif name == "FORCE":
                records = []
                self.assertEqual(float(fields(lines[:1])[2]), 50.0)
                for line in lines[1:]:
                    if line[8:16].strip():
                        records.append([])
                    records[-1].extend(map(int, fields([line])))
                self.assertEqual(len(records), int(fields(lines[:1])[-1]))
                for record in records:
                    self.assertEqual(len(record), size + 1)
                    for local, value in enumerate(record[1:]):
                        if value:
                            decoded_loads.add(location(record[0], local))
        expected = np.zeros(shape, dtype=np.uint8)
        expected[tuple(slice(0, n) for n in occupancy.shape)] = occupancy
        np.testing.assert_array_equal(rebuilt, expected)
        self.assertEqual(decoded_loads, loads)
        for name, lines in parsed:
            if name == "COARSE":
                eid, mtype, _ = map(int, fields(lines))
                block = [rebuilt[location(eid, local)] for local in range(size)]
                self.assertEqual(mtype, 0 if not any(block) else 1 if all(block) else 2)
        spc = list(map(int, fields(next(lines for name, lines in parsed if name == "SPC1"))))
        expected_nodes = {1 + x + (coarse[0] + 1) * y + (coarse[0] + 1) * (coarse[1] + 1) * z
                          for x in (0, 1) for y in (0, 1) for z in (0, 1)}
        self.assertEqual(set(spc[2:]), expected_nodes)

    def test_writers_preserve_voxels_spacing_loads_and_constraints(self):
        for fine in ((5, 5, 5), (10, 10, 10), (2, 3, 4)):
            for stream in (False, True):
                with self.subTest(fine=fine, stream=stream):
                    args = (self.occ, (4., 5., 6.), (.25, .5, .75), {(0, 0, 0)},
                            self.loads, fine, "123", 1, 50., (0., 0., -1.), 1, "test")
                    if stream:
                        path = self.root / "test.vgm"
                        vgm.write_vgm_stream(path, *args)
                        text = path.read_text()
                    else:
                        text = "\n".join(vgm.build_vgm(*args))
                    self.check_geometry_and_loads(text, fine)

    def test_presets_defaults_and_invalid_values(self):
        parser = argparse.ArgumentParser()
        vgm.add_fine_arguments(parser)
        for argv, expected in (([], (5, 5, 5)), (["--fine-preset", "101010"], (10, 10, 10)),
                               (["--fine-preset", "555"], (5, 5, 5)),
                               (["--fine", "10", "10", "10"], (10, 10, 10))):
            self.assertEqual(vgm.resolve_fine(parser.parse_args(argv), parser), expected)
        for argv in (["--fine", "0", "10", "10"], ["--fine", "-1", "5", "5"],
                     ["--fine-preset", "101010", "--fine", "5", "5", "5"]):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                vgm.resolve_fine(parser.parse_args(argv), parser)

    def write_occupancy(self):
        path = self.root / "occupancy.vtk"
        path.write_text("# vtk DataFile Version 3.0\ntest\nASCII\nDATASET STRUCTURED_POINTS\n"
                        "DIMENSIONS 22 12 13\nORIGIN 4 5 6\nSPACING 0.25 0.5 0.75\n"
                        f"CELL_DATA {self.occ.size}\nSCALARS occupancy unsigned_char\nLOOKUP_TABLE default\n"
                        + " ".join(map(str, self.occ.ravel(order="F"))) + "\n")
        return path

    def write_selection(self, name, voxels):
        path = self.root / name
        points = list(voxels)
        text = ("# vtk DataFile Version 3.0\nselection\nASCII\nDATASET POLYDATA\n"
                f"CELL_DATA {len(points)}\n")
        for axis, label in enumerate(("voxel_i", "voxel_j", "voxel_k")):
            text += f"SCALARS {label} int\nLOOKUP_TABLE default\n"
            text += " ".join(str(point[axis]) for point in points) + "\n"
        path.write_text(text)
        return path

    def test_vtk_cli_101010_and_explicit_fine(self):
        occ = self.write_occupancy()
        constraint = self.write_selection("constraint.vtk", {(0, 0, 0)})
        load = self.write_selection("load.vtk", self.loads)
        outputs = []
        for fine_args in (["--fine-preset", "101010"], ["--fine", "10", "10", "10"]):
            output = self.root / "cli.vgm"
            result = subprocess.run([sys.executable, str(TOOLKIT / "vtk_selection_to_vgm.py"),
                                     "--occupancy", str(occ), "--constraint-selection", str(constraint),
                                     "--load-selection", str(load), "--output", str(output),
                                     "--force", "50", *fine_args], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(output.read_text())
            self.check_geometry_and_loads(outputs[-1], (10, 10, 10))
        self.assertEqual(*outputs)

    def test_stl_cli_101010(self):
        occ = self.write_occupancy()
        # Small triangles around the centers of two occupied fine cells.
        paths = []
        for name, center in (("constraint.stl", (4.125, 5.25, 6.375)),
                             ("load.stl", (9.125, 10.25, 14.625))):
            x, y, z = center
            path = self.root / name
            path.write_text(f"solid face\nfacet normal 0 0 1\nouter loop\n"
                            f"vertex {x-.01} {y-.01} {z}\nvertex {x+.01} {y-.01} {z}\n"
                            f"vertex {x} {y+.01} {z}\nendloop\nendfacet\nendsolid face\n")
            paths.append(path)
        output = self.root / "surface.vgm"
        result = subprocess.run([sys.executable, str(TOOLKIT / "stl_surface_to_vgm.py"),
                                 "--occupancy", str(occ), "--constraint-stl", str(paths[0]),
                                 "--load-stl", str(paths[1]), "--output", str(output),
                                 "--fine-preset", "101010", "--force", "50", "--tolerance", "0.01",
                                 "--workers", "1"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.check_geometry_and_loads(output.read_text(), (10, 10, 10), loads={(20, 10, 11)})

    def test_gui_preset_and_custom_fields_stay_in_sync(self):
        from gui import VoxelGui
        # A hidden Tk window exercises the actual widgets without starting a user session.
        app = VoxelGui()
        app.withdraw()
        self.addCleanup(app.destroy)
        self.assertEqual(app.vgm_fine_preset_var.get(), "5x5x5")
        app.vgm_fine_preset_var.set("10x10x10")
        app._apply_fine_preset()
        self.assertEqual([v.get() for v in (app.vgm_fine_x_var, app.vgm_fine_y_var,
                                          app.vgm_fine_z_var)], ["10"] * 3)
        app.vgm_occ_var.set(str(self.write_occupancy()))
        app.vgm_constraint_var.set(str(self.write_selection("constraint.vtk", {(0, 0, 0)})))
        app.vgm_load_var.set(str(self.write_selection("load.vtk", self.loads)))
        app.vgm_output_var.set(str(self.root / "gui.vgm"))
        self.assertEqual(app._parse_vgm_config()["fine"], (10, 10, 10))
        app.vgm_fine_x_var.set("7")
        self.assertEqual(app.vgm_fine_preset_var.get(), "自定义")
        app.vgm_fine_preset_var.set("5x5x5")
        app._apply_fine_preset()
        self.assertEqual(app.vgm_fine_preset_var.get(), "5x5x5")


if __name__ == "__main__":
    unittest.main()
