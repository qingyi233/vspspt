# STP Voxel VTK Toolkit

这个文件夹是整理后的可用版本，包含 GUI、命令行工具、时间对比脚本、C++ 源码、已编译 C++ 运行库。输入支持 `.stp/.step/.stl`。

## 目录

- `gui.py`：图形界面，选择 `.stp/.step/.stl` 后自动输出体素 VTK。
- `stp_to_voxel_vtk.py`：命令行转换工具，适合批处理和亿级体素。
- `benchmark_time_comparison.py`：C++ 与 Python 体素化耗时对比。
- `pyvoxel_runtime/`：最小运行时，包含 `cvoxelise.pyd`。
- `cpp_cvoxelise/`：C++/pybind11 源码和 Visual Studio 工程。
- `third_party/eigen/`：C++ 编译需要的 Eigen 头文件。
- `docs/TIME_COMPARISON.md`：已跑过的时间对比报告。
- `output/`：默认输出目录。

## 先安装依赖

```powershell
cd "D:\File\Project\Section\体素化代码\stp_voxel_toolkit"
python -m pip install numpy gmsh joblib pybind11
```

## GUI 使用

```powershell
python gui.py
```

建议设置：

- `Target cells`：`100000000` 或 `1e8`，生成约 1 亿总体素单元。
- 勾选 `Use C++ voxelizer`。
- 勾选 `Write surface VTK`，这个文件能在 ParaView 中直接看到结构。
- 勾选 `Reuse STL`，第二次运行会跳过 STEP 转 STL，速度明显更快。

GUI 会显示各模块耗时、占比和累计时间。

## 命令行使用

快速预览：

```powershell
python stp_to_voxel_vtk.py --input "..\test1.stp" --resolution 400 --surface --reuse-stl
```

亿级总体素：

```powershell
python stp_to_voxel_vtk.py --input "..\test1.stp" --target-cells 100000000 --surface --reuse-stl
```

也可以用科学计数法：

```powershell
python stp_to_voxel_vtk.py --input "..\test1.stp" --target-cells 1e8 --surface --reuse-stl
```

输出文件：

- `*_occupancy_*.vtk`：结构化体素占据场，适合后续计算；ParaView 中需要 Threshold `occupancy = 1`。
- `*_surface_*.vtk`：从占据体素提取出的外表面，适合直接显示。

## STEP 三角化失败

如果看到类似 `Impossible to mesh periodic surface` 的错误，说明失败发生在 `STEP -> STL` 三角化阶段，通常是 CAD 模型存在周期面、坏面、退化边或缝合问题，不是体素化或 C++ 代码的问题。

处理办法：

- 在 CAD 软件中修复/缝合几何后重新导出 STEP。
- 直接从 CAD 软件导出 STL，然后在本工具中选择 STL 文件继续体素化。
- 对复杂装配先拆分或简化，再转换。

## 时间对比

```powershell
python benchmark_time_comparison.py --input "..\test1.stp" --resolution 400
```

已测结果见：

```text
docs\TIME_COMPARISON.md
```

当前机器上，`resolution=400` 时：

- C++ 体素化：`1.140 s`
- Python 体素化：`11.682 s`
- C++ 加速约 `10.25x`

## C++ 说明

默认运行已经使用 `pyvoxel_runtime\cvoxelise.pyd`，不用重新编译。

需要重新编译时看：

```text
cpp_cvoxelise\BUILD.md
```

注意：这个环境中 C++ 的 `parallel=True` 不稳定，工具默认使用 C++ 串行路径 `parallel=False`，速度仍明显快于 Python。
