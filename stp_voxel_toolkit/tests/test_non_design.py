"""Check that non-design geometry selects only original analysis cells."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from test_model_generation import cube_mesh
from test_vgm_subdivisions import cards, fields
from non_design_to_vgm import select_non_design
from stp_to_voxel_vtk import write_ascii_stl


class NonDesignTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.analysis = self.root / "analysis.vgm"
        self.geometry = self.root / "region.stl"
        self.output = self.root / "mask.vgm"
        self.origin = (10., 20., 30.)
        self.box = "BOX          6.0     8.0     1.0       3       2       1\n               2       2       2"
        self.source = "\n".join((
            "$ ORIGIN 10 20 30", self.box,
            "MATB           4       0       2       0       1       1       0       2+M4",
            "+M4                    0",
            "COARSE         1       1       1",
            "COARSE         2       2       1",
            "COARSE         3       0       1",
            "COARSE         4       2       1",
            "COARSE         5       1       1",
            "COARSE         6       0       1",
            "MATB           2       1       0       2       1       0       2       0",
            "                       1",
            "ENDDATA\n",
        ))
        self.analysis.write_text(self.source, encoding="ascii")
        self.write_region((8.2, 20.3, 29.8), (14.7, 25.2, 31.2))

    def write_region(self, low, high, closed=True):
        unit_cube = (cube_mesh().transpose(0, 2, 1) - (2, 3, 4)) / 2
        triangles = low + unit_cube * np.subtract(high, low)
        if not closed:
            triangles = triangles[:-1]
        write_ascii_stl(self.geometry, "region", triangles.reshape(-1, 3).tolist())

    def test_selects_original_eids_and_fine_indices_in_translated_box(self):
        stats = select_non_design(self.analysis, self.geometry, self.output)
        self.assertEqual(stats, dict(candidate_cells=42, selected_cells=19, excluded_cells=23,
                                     selected_coarse_cells=4))
        exported = self.output.read_text(encoding="ascii")
        self.assertIn("$ ORIGIN 10 20 30\n", exported)
        self.assertIn(self.box, exported)
        expected = {1: list(range(8)), 2: [0, 2, 3, 5, 7], 4: [1, 4], 5: [0, 1, 4, 5]}
        coarse = {}
        selected = {}
        for name, lines in cards(exported):
            if name == "COARSE":
                eid, cell_type, optimizable = map(int, fields(lines))
                self.assertEqual(optimizable, 0)
                coarse[eid] = cell_type
            elif name == "MATB":
                eid, *mask = map(int, fields(lines))
                self.assertEqual(len(mask), 8)
                selected[eid] = [i for i, value in enumerate(mask) if value]
        self.assertEqual(coarse, {1: 1, 2: 2, 3: 0, 4: 2, 5: 2, 6: 0})
        self.assertEqual(selected, expected)
        self.assertEqual(self.analysis.read_text(encoding="ascii"), self.source)

    def test_gui_accepts_explicit_origin_and_defaults_output_path(self):
        from gui import VoxelGui

        self.analysis.write_text(self.source.partition("\n")[2], encoding="ascii")
        app = VoxelGui()
        app.withdraw()
        self.addCleanup(app.destroy)
        app.notebook.select(2)
        app._tab_changed()
        self.assertEqual(app.non_design_button.winfo_manager(), "pack")
        app.non_design_analysis_var.set(str(self.analysis))
        app.non_design_geometry_var.set(str(self.geometry))
        app.non_design_origin_var.set("10,20,30")
        with patch.object(app, "_begin_task") as begin:
            app.non_design_button.invoke()
        target, config, _status = begin.call_args.args
        target(config)
        app._poll_queue()
        self.assertTrue(self.analysis.with_name("analysis_NonDesignDomain_subset.vgm").is_file())
        self.assertEqual(app.status_var.get(), "非设计域已提取")
        self.assertIn("选中 19 个实体细单元", app.log_text.get("1.0", "end"))

    def test_invalid_geometry_and_output_leave_inputs_intact(self):
        with self.assertRaisesRegex(ValueError, "differs from"):
            select_non_design(self.analysis, self.geometry, self.output, origin=(0., 0., 0.))
        with self.assertRaisesRegex(ValueError, "must be separate"):
            select_non_design(self.analysis, self.geometry, self.analysis)
        self.write_region((14.1, 20.3, 29.8), (15.8, 23.8, 31.2))
        with self.assertRaisesRegex(ValueError, "did not select any solid"):
            select_non_design(self.analysis, self.geometry, self.output)
        self.write_region((8.2, 20.3, 29.8), (14.7, 25.2, 31.2), closed=False)
        with self.assertRaisesRegex(ValueError, "closed triangle mesh"):
            select_non_design(self.analysis, self.geometry, self.output)
        self.assertFalse(self.output.exists())
        self.assertEqual(self.analysis.read_text(encoding="ascii"), self.source)
