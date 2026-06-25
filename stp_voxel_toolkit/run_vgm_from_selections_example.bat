@echo off
cd /d "%~dp0"
python vtk_selection_to_vgm.py ^
  --occupancy "output\your_occupancy.vtk" ^
  --constraint-selection "constraint_selection.vtk" ^
  --load-selection "load_selection.vtk" ^
  --output "output\model.vgm" ^
  --fine 5 5 5 ^
  --spc-dof 123 ^
  --force 1000 ^
  --direction 0,0,-1
