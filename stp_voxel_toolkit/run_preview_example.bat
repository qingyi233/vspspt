@echo off
cd /d "%~dp0"
python stp_to_voxel_vtk.py --input "..\test1.stp" --resolution 400 --surface --reuse-stl
pause
