"""End-to-end checks for model-only VGM, custom boxes and mesh statistics."""

import contextlib
import io
import itertools
import json
from pathlib import Path
import queue
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

TOOLKIT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLKIT))
import vtk_selection_to_vgm as vgm
from vgm_stats import active_node_count, mesh_statistics
from voxel_grid import resolve_grid, validate_box


def cube_mesh():
    points = np.array([(2, 3, 4), (4, 3, 4), (4, 5, 4), (2, 5, 4),
                       (2, 3, 6), (4, 3, 6), (4, 5, 6), (2, 5, 6)], float)
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    triangles = [points[list(tri)] for a, b, c, d in quads for tri in ((a, b, c), (a, c, d))]
    return np.array(triangles).transpose(0, 2, 1)


def reference_nodes(occupancy):
    return len({tuple(int(i) + d for i, d in zip(cell, delta))
                for cell in np.argwhere(occupancy)
                for delta in itertools.product((0, 1), repeat=3)})


class GenerationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def occupancy_file(self):
        from stp_to_voxel_vtk import write_occupancy_vtk
        path = self.root / "occupancy.vtk"
        write_occupancy_vtk(path, np.ones((2, 2, 2), np.uint8), (2., 3., 4.), (1., 1., 1.))
        return path

    def cube_file(self):
        from stp_to_voxel_vtk import write_ascii_stl
        path = self.root / "cube.stl"
        write_ascii_stl(path, "cube", cube_mesh().transpose(0, 2, 1).reshape(-1, 3).tolist())
        return path

    def test_tangent_rays_at_grid_edges_preserve_solid_cells(self):
        from gui import pvoxel, voxel

        vertices = np.array([(0, 0, 0), (0, 2, 0), (0, 0, 2), (2, 0, 0)], dtype=float)
        triangles = vertices[[(0, 1, 2), (0, 2, 3), (0, 3, 1), (1, 3, 2)]].transpose(0, 2, 1)
        bounds = np.array([[0., 0., 0.], [2., 2., 2.]])
        # The Z rays at x + y = 2 hit the tetrahedron at just one point.
        centers = 0.25 + 0.5 * np.indices((4, 4, 4))
        expected = centers.sum(axis=0) < 2.
        for backend in (pvoxel, voxel):
            with self.subTest(backend=backend.__name__), contextlib.redirect_stdout(io.StringIO()):
                grid = backend.Voxel(4, 4, 4).VOXELISE(triangles, ray="z", parallel=False, bounds=bounds)
                if grid.ndim == 4:
                    grid = grid[..., 0]
                np.testing.assert_array_equal(grid, expected)

    def run_count_button(self, app, button):
        with patch.object(app, "_begin_task") as begin, patch("gui.messagebox.showerror") as error:
            button.invoke()
            error.assert_not_called()
            worker, request, _ = begin.call_args.args
        with contextlib.redirect_stdout(io.StringIO()), patch("gui.messagebox.showerror") as error:
            worker(request)
            app._poll_queue()
            error.assert_not_called()

    def test_node_counts_match_independent_corner_union(self):
        rng = np.random.default_rng(12)
        for shape in ((1, 1, 1), (2, 3, 4), (4, 1, 2), (3, 5, 2)):
            for occupancy in (np.zeros(shape, np.uint8), np.ones(shape, np.uint8),
                              rng.integers(0, 2, size=shape, dtype=np.uint8)):
                with self.subTest(shape=shape, occupied=int(occupancy.sum())):
                    self.assertEqual(active_node_count(occupancy), reference_nodes(occupancy))
        occupancy = np.zeros((4, 4, 4), np.uint8)
        occupancy[0, 0, 0] = occupancy[1, 1, 1] = 1  # Two cells sharing just one node.
        stats = mesh_statistics(occupancy, (2, 2, 2), (0, 0, 0), (1, 1, 1))
        self.assertEqual(stats["solid_fine_nodes"], 15)
        self.assertEqual(stats["solid_fine_dofs"], 45)
        self.assertEqual(stats["solid_fine_elements"], 2)
        self.assertEqual(stats["active_coarse_elements"], 1)
        self.assertEqual(stats["coarse_grid_elements"], 8)
        self.assertEqual(stats["active_coarse_nodes"], 8)

    def test_optional_boundaries_have_no_dangling_cards_or_references(self):
        for constraints, loads, stream in itertools.product((False, True), repeat=3):
            with self.subTest(constraints=constraints, loads=loads, stream=stream):
                occupancy = np.ones((2, 2, 2), np.uint8)
                occupancy[0, 0, 1] = 0
                candidates = {(0, 0, 1), (1, 1, 1)} if loads else {(0, 0, 1)}
                args = (occupancy, (0., 0., 0.), (1., 1., 1.),
                        {(0, 0, 0)} if constraints else set(), candidates,
                        (2, 2, 2), "13", 1, 50., (0., 0., -1.), 1, "test")
                if stream:
                    path = self.root / "model.vgm"
                    stats = vgm.write_vgm_stream(path, *args)
                    text = path.read_text()
                    self.assertEqual(json.loads(path.with_suffix(".stats.json").read_text()), stats)
                    self.assertEqual(stats["constrained_coarse_dofs"], 16 if constraints else 0)
                else:
                    text = "\n".join(vgm.build_vgm(*args))
                names = {line[:8].strip() for line in text.splitlines()}
                for name in ("SPC1", "SPCADD"):
                    self.assertEqual(name in names, constraints)
                for name in ("FORCE", "LOAD"):
                    self.assertEqual(name in names, loads)
                self.assertEqual("   SPC = 2" in text, constraints)
                self.assertEqual("   LOAD = 21" in text, loads)
                self.assertIn("MATB", names)
                self.assertTrue(text.endswith("ENDDATA\n") or text.endswith("ENDDATA"))

    def test_vtk_and_stl_clis_accept_no_boundaries(self):
        occ = self.occupancy_file()
        for script in ("vtk_selection_to_vgm.py", "stl_surface_to_vgm.py"):
            with self.subTest(script=script):
                path = self.root / (script + ".vgm")
                result = subprocess.run([sys.executable, str(TOOLKIT / script), "--occupancy", str(occ),
                                         "--output", str(path)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                stats = json.loads(path.with_suffix(".stats.json").read_text())
                self.assertEqual(stats["solid_fine_elements"], 8)
                self.assertEqual(stats["solid_fine_nodes"], 27)
                self.assertEqual(stats["solid_fine_dofs"], 81)
                self.assertIn("81", result.stdout)

    def test_box_rejects_every_uncontained_face_and_invalid_values(self):
        mesh = cube_mesh()
        for axis in range(3):
            for minimum_face in (False, True):
                origin = np.array([2., 3., 4.])
                size = np.array([2., 2., 2.])
                if minimum_face:
                    origin[axis] += 0.01
                else:
                    size[axis] -= 0.01
                with self.subTest(axis=axis, minimum_face=minimum_face), self.assertRaisesRegex(ValueError, "未完全包含"):
                    resolve_grid(mesh, 64, None, origin, size)
        for origin, size in (((0, 0, 0), (0, 2, 2)), ((0, 0, 0), (-1, 2, 2)),
                             ((float("nan"), 0, 0), (4, 4, 4)), ((0, 0, 0), (float("inf"), 2, 2))):
            with self.subTest(origin=origin, size=size), self.assertRaises(ValueError):
                validate_box(origin, size)
        with self.assertRaisesRegex(ValueError, "同时指定"):
            resolve_grid(mesh, 64, None, (1, 2, 3), None)

    def test_box_sampling_matches_world_coordinates_in_both_backends(self):
        from stp_to_voxel_vtk import pvoxel, voxel
        expected = np.zeros((4, 4, 4), np.uint8)
        expected[1:3, 1:3, 1:3] = 1
        bounds = np.array([[1., 2., 3.], [5., 6., 7.]])
        for backend, ray in itertools.product((pvoxel, voxel), ("xyz", "z")):
            with self.subTest(backend=backend.__name__, ray=ray), contextlib.redirect_stdout(io.StringIO()):
                obj = backend.Voxel(4, 4, 4)
                grid = np.asarray(obj.VOXELISE(cube_mesh(), ray=ray, parallel=False, bounds=bounds))
                if grid.ndim == 4:
                    grid = grid[..., 0]
                np.testing.assert_array_equal(grid, expected)
                for axis, values in enumerate((obj.gridCOx, obj.gridCOy, obj.gridCOz)):
                    np.testing.assert_allclose(values, bounds[0, axis] + np.arange(4) + 0.5)

    def test_direct_cli_vgm_keeps_exact_custom_box_and_rejects_clipping(self):
        stl = self.cube_file()
        for pure_python in (False, True):
            with self.subTest(pure_python=pure_python):
                path = self.root / f"cube_{pure_python}.vgm"
                argv = [sys.executable, str(TOOLKIT / "stp_to_voxel_vtk.py"), "--input", str(stl),
                        "--output", str(self.root / "outputs"), "--resolution", "4", "--vgm", str(path),
                        "--no-occupancy", "--fine", "2", "2", "2", "--box-origin", "1", "2", "3",
                        "--box-size", "4", "4", "4"]
                result = subprocess.run(argv + (["--python"] if pure_python else []), capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                stats = json.loads(path.with_suffix(".stats.json").read_text())
                self.assertEqual(stats["origin"], [1, 2, 3])
                self.assertEqual(stats["box_size"], [4, 4, 4])
                self.assertEqual(stats["solid_fine_elements"], 8)
                self.assertEqual(stats["solid_fine_nodes"], 27)
                path.unlink()
                argv[-3:] = ["2", "2", "2"]  # The box no longer contains the cube.
                result = subprocess.run(argv + (["--python"] if pure_python else []), capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("未完全包含", result.stderr)
                self.assertFalse(path.exists())
        self.assertFalse(list((self.root / "outputs").glob("*.vtk")))
        dims, origin, spacing = resolve_grid(cube_mesh(), 125, None, (1, 2, 3), (4, 4, 4), fine=(3, 3, 3))
        np.testing.assert_array_equal(np.array(dims) % 3, 0)
        np.testing.assert_allclose(np.array(dims) * spacing, (4, 4, 4))

    def test_gui_model_only_pipeline_and_optional_selection_fields(self):
        from gui import VoxelGui, run_pipeline
        app = VoxelGui()
        app.withdraw()
        self.addCleanup(app.destroy)
        app.step_var.set(str(self.occupancy_file()))
        app.out_var.set(str(self.root))
        app.vgm_force_var.set("")
        app.vgm_direction_var.set("")
        app.vgm_spc_dof_var.set("")
        config = app._parse_vgm_config()
        self.assertIsNone(config["constraint"])
        self.assertIsNone(config["load"])
        app._run_vgm_worker(config)
        messages = list(app.queue.queue)
        self.assertFalse(any(item[0] == "error" for item in messages), messages)
        self.assertTrue(any(item[0] == "vgm_statistics" for item in messages))
        with patch("gui.messagebox.showinfo") as success_dialog:
            app._poll_queue()
            success_dialog.assert_not_called()
        self.assertIn("81", app.vgm_stats_var.get())
        self.assertEqual(app.status_var.get(), "生成完成")

        stl = self.cube_file()
        app.step_var.set(str(stl))
        app.write_occ_var.set(False)
        app.write_surface_var.set(False)
        app.write_vgm_var.set(True)
        app.vgm_fine_x_var.set("2")
        app.vgm_fine_y_var.set("2")
        app.vgm_fine_z_var.set("2")
        app.custom_box_var.set(True)
        app.box_origin_var.set("1,2,3")
        app.box_size_var.set("4,4,4")
        app.resolution_var.set("4")
        app.grid_mode_var.set("最长轴格数")
        app.out_var.set(str(self.root / "gui_output"))
        config = app._parse_config()
        self.assertTrue(config.write_vgm)
        q = queue.Queue()
        with contextlib.redirect_stdout(io.StringIO()):
            run_pipeline(config, q)
        paths = list(config.out_dir.glob("*.vgm"))
        self.assertEqual(len(paths), 1)
        stats = json.loads(paths[0].with_suffix(".stats.json").read_text())
        self.assertEqual(stats["solid_fine_dofs"], 81)
        self.assertTrue(any(item[0] == "done" and item[1] for item in q.queue))
        config.box_size = (2., 2., 2.)
        with self.assertRaisesRegex(ValueError, "未完全包含"):
            run_pipeline(config, queue.Queue())

    def test_step_to_model_only_vgm_end_to_end(self):
        import gmsh
        step = self.root / "cube.step"
        gmsh.initialize(interruptible=False)
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("cube")
            gmsh.model.occ.addBox(2, 3, 4, 2, 2, 2)
            gmsh.model.occ.synchronize()
            gmsh.write(str(step))
        finally:
            gmsh.finalize()
        path = self.root / "step_result.vgm"
        result = subprocess.run([
            sys.executable, str(TOOLKIT / "stp_to_voxel_vtk.py"), "--input", str(step),
            "--output", str(self.root / "step_output"), "--resolution", "4", "--mesh-size", "1",
            "--vgm", str(path), "--no-occupancy", "--box-origin", "1", "2", "3",
            "--box-size", "4", "4", "4", "--fine", "2", "2", "2"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        stats = json.loads(path.with_suffix(".stats.json").read_text())
        self.assertEqual(stats["solid_fine_elements"], 8)
        self.assertEqual(stats["solid_fine_nodes"], 27)
        self.assertEqual(stats["solid_fine_dofs"], 81)

    def test_gui_unified_entry_routes_cad_and_occupancy_and_ignores_inactive_settings(self):
        from gui import VoxelGui
        app = VoxelGui()
        app.withdraw()
        self.addCleanup(app.destroy)
        self.assertEqual([app.notebook.tab(tab, "text") for tab in app.notebook.tabs()],
                         ["模型生成", "有限元求解"])
        app.out_var.set(str(self.root))
        app.step_var.set(str(self.cube_file()))
        app.grid_mode_var.set("最长轴格数")
        app.resolution_var.set("4")
        app.target_var.set("ignored")
        app.box_origin_var.set("ignored")  # The automatic box is selected.
        config = app._parse_config()
        self.assertEqual(config.longest_resolution, 4)
        self.assertIsNone(config.box_origin)
        with patch.object(app, "_begin_task") as begin:
            app._start()
            self.assertEqual(begin.call_args.args[0].__name__, "_run_worker")
        app.grid_mode_var.set("目标体素数")
        app.target_var.set("64")
        app.resolution_var.set("ignored")
        self.assertIsNone(app._parse_config().longest_resolution)

        app.step_var.set(str(self.occupancy_file()))
        self.assertEqual(app.grid_frame.winfo_manager(), "")
        self.assertEqual(app.outputs_frame.winfo_manager(), "")
        self.assertEqual(app.vgm_frame.winfo_manager(), "pack")
        self.assertEqual(app.constraint_fields.winfo_manager(), "")
        self.assertEqual(app.load_fields.winfo_manager(), "")
        app.vgm_workers_var.set("ignored")
        app.vgm_spc_dof_var.set("ignored")
        app.vgm_direction_var.set("ignored")
        app.vgm_force_var.set("ignored")
        with patch.object(app, "_begin_task") as begin:
            app._start()
            self.assertEqual(begin.call_args.args[0].__name__, "_run_vgm_worker")
            self.assertEqual(begin.call_args.args[1]["output"], self.root / "occupancy.vgm")
        app.vgm_fine_preset_var.set("自定义")
        self.assertEqual(app.fine_fields.winfo_manager(), "pack")
        app.vgm_fine_preset_var.set("10x10x10")
        app._apply_fine_preset()
        self.assertEqual(app.fine_fields.winfo_manager(), "")

        app.step_var.set(str(self.cube_file()))
        app.write_vgm_var.set(False)
        app.write_occ_var.set(True)
        self.assertEqual(app.vgm_frame.winfo_manager(), "")
        self.assertEqual(app.grid_frame.winfo_manager(), "pack")
        self.assertIsNone(app._parse_config().vgm_options)

    def test_gui_direct_cad_generation_applies_optional_boundaries(self):
        from gui import VoxelGui, run_pipeline
        app = VoxelGui()
        app.withdraw()
        self.addCleanup(app.destroy)
        app.step_var.set(str(self.cube_file()))
        app.out_var.set(str(self.root / "with_boundaries"))
        app.grid_mode_var.set("最长轴格数")
        app.resolution_var.set("4")
        app.custom_box_var.set(True)
        app.box_origin_var.set("1,2,3")
        app.box_size_var.set("4,4,4")
        for var in (app.vgm_fine_x_var, app.vgm_fine_y_var, app.vgm_fine_z_var):
            var.set("2")
        selection = self.root / "selection.vtk"
        selection.write_text("# vtk DataFile Version 3.0\nselection\nASCII\nDATASET POLYDATA\nCELL_DATA 1\n"
                             + "".join(f"SCALARS voxel_{axis} int\nLOOKUP_TABLE default\n1\n" for axis in "ijk"))
        app.vgm_constraint_var.set(str(selection))
        app.vgm_load_var.set(str(selection))
        self.assertEqual(app.constraint_fields.winfo_manager(), "pack")
        self.assertEqual(app.load_fields.winfo_manager(), "pack")
        self.run_count_button(app, app.load_count_button)
        self.assertEqual(app.grid_count_var.get(), "实体体素 8")
        self.assertEqual(app.load_count_var.get(), "预计命中 1 个实体细单元")
        self.assertFalse((self.root / "with_boundaries").exists())
        with contextlib.redirect_stdout(io.StringIO()):
            run_pipeline(app._parse_config(), queue.Queue())
        text = (self.root / "with_boundaries" / "cube.vgm").read_text()
        self.assertIn("   SPC = 2", text)
        self.assertIn("   LOAD = 21", text)
        self.assertIn("FORCE", text)

    def test_gui_counts_solid_load_cells_and_exports_requested_total(self):
        from gui import FORCE_TOTAL, VoxelGui
        from stp_to_voxel_vtk import write_occupancy_vtk
        from test_vgm_subdivisions import cards, fields

        app = VoxelGui()
        app.withdraw()
        self.addCleanup(app.destroy)
        occupancy = np.zeros((3, 2, 2), np.uint8)
        occupancy[0, 0, 0] = occupancy[1, 0, 0] = occupancy[2, 1, 1] = 1
        path = self.root / "sparse.vtk"
        write_occupancy_vtk(path, occupancy, (0., 0., 0.), (1., 1., 1.))
        app.step_var.set(str(path))
        app.out_var.set("")
        self.run_count_button(app, app.grid_count_button)
        self.assertEqual(app.grid_count_var.get(), "实体体素 3")

        selection = self.root / "load.vtk"
        header = "# vtk DataFile Version 3.0\nload\nASCII\nDATASET POLYDATA\n"
        scalar = "SCALARS voxel_eid int\nLOOKUP_TABLE default\n"
        selection.write_text(header + "CELL_DATA 5\n" + scalar + "1 2 1 3 13\n")
        app.vgm_load_var.set(str(selection))
        app.vgm_force_var.set("")  # Counting does not require load magnitude.
        self.run_count_button(app, app.load_count_button)
        self.assertEqual(app.load_count_var.get(), "预计命中 2 个实体细单元")
        self.assertFalse(list(self.root.glob("*.vgm")))

        app.out_var.set(str(self.root))
        self.assertEqual(app.load_count_var.get(), "预计命中 2 个实体细单元")
        app.vgm_force_mode_var.set(FORCE_TOTAL)
        app.vgm_force_var.set("90")
        app.vgm_direction_var.set("0,0,-2")
        app._run_vgm_worker(app._parse_vgm_config())
        output = self.root / "sparse.vgm"
        force = next(lines for name, lines in cards(output.read_text()) if name == "FORCE")
        self.assertEqual(float(fields(force[:1])[2]), 45.)
        self.assertEqual([float(value) for value in fields(force[:1])[3:6]], [0., 0., -1.])
        self.assertEqual(list(map(int, fields(force[1:])))[1:], [1, 1] + [0] * 123)

        empty_selection = self.root / "void_load.vtk"
        empty_selection.write_text(header + "CELL_DATA 1\n" + scalar + "3\n")
        app.vgm_load_var.set(str(empty_selection))
        app.out_var.set(str(self.root / "empty"))
        app._run_vgm_worker(app._parse_vgm_config())
        self.assertTrue(any(item[0] == "error" and "未命中实体细单元" in item[1] for item in app.queue.queue))
        self.assertFalse((self.root / "empty" / "sparse.vgm").exists())


if __name__ == "__main__":
    unittest.main()
