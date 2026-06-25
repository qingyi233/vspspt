# STP Voxel VTK Toolkit

这个文件夹是整理后的可用版本，包含 GUI、命令行工具、时间对比脚本、C++ 源码、已编译 C++ 运行库。输入支持 `.stp/.step/.stl`。

## 目录

- `gui.py`：图形界面，选择 `.stp/.step/.stl` 后自动输出体素 VTK。
- `stp_to_voxel_vtk.py`：命令行转换工具，适合批处理和亿级体素。
- `vtk_selection_to_vgm.py`：读取体素 VTK 和 ParaView 保存的两个选区 VTK，生成自定义 VGM。
- `stl_surface_to_vgm.py`：读取体素 VTK 和 UG/ParaView 导出的两个曲面 STL，生成自定义 VGM。
- `run_dlutfem.py`：调用 `DLUTFEM.exe` 计算 VGM，并自动绕开中文路径问题。
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

如果使用源码版压缩包，解压后先安装上面的 Python 依赖，再双击 `run_stp_voxel_toolkit.bat`。

如果使用便携版压缩包 `stp_voxel_toolkit_portable.zip`，解压后直接双击 `run_stp_voxel_toolkit.bat`。便携版自带 Python、numpy、gmsh、joblib、体素化运行库和 `solver_bin\DLUTFEM.exe`，不需要在对方电脑安装 Python。

如果给外部用户分发，优先使用安装器 `STPVoxelToolkit_Setup.exe`：双击后会安装到当前用户的 `%LOCALAPPDATA%\STPVoxelToolkit`，并创建桌面和开始菜单快捷方式，不需要管理员权限。

## GUI 使用

```powershell
python gui.py
```

建议设置：

- `Target cells`：`100000000` 或 `1e8`，生成约 1 亿总体素单元。
- 勾选 `Use C++ voxelizer`。
- 勾选 `Write surface VTK`，这个文件能在 ParaView 中直接看到结构。
- 勾选 `Reuse STL`，第二次运行会跳过 STEP 转 STL，速度明显更快。

GUI 现在有三个分页：

- `体素化生成 VTK`：把 STP/STEP/STL 转成体素 VTK。
- `选区生成 VGM`：用 occupancy VTK 加两个选区 VTK/STL 生成 VGM。
- `VGM 有限元求解`：调用 `DLUTFEM.exe` 计算 VGM，并把结果 `.h5/.vtk/.db/.log` 回拷到指定目录。

体素化页会显示各模块耗时、占比和累计时间；VGM 生成和有限元求解会在日志窗口持续输出状态。

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
- `*_surface_*.vtk`：从占据体素提取出的外表面，适合直接显示。该文件的 cell data 带有 `voxel_i`、`voxel_j`、`voxel_k`、`voxel_eid`、`face_id`，用于 ParaView 选区后回写 VGM。

## ParaView 选区生成 VGM

第一版推荐在 ParaView 中选择 `*_surface_*.vtk`，不要直接选择 STL 或 STP。表面 VTK 的每个面都带有所属体素编号，后续能稳定映射到 VGM。

1. 重新生成带映射字段的体素文件：

```powershell
python stp_to_voxel_vtk.py --input "..\test2_stp.stl" --target-cells 100000000 --surface
```

2. 在 ParaView 打开生成的 `*_surface_*.vtk`。
3. 选择约束区域，使用 `Extract Selection` 提取选区，然后 `File -> Save Data` 保存为 legacy `.vtk`，例如：

```text
constraint_selection.vtk
```

4. 选择载荷区域，同样 `Extract Selection` 后保存为：

```text
load_selection.vtk
```

5. 用原始 occupancy VTK 和两个选区 VTK 生成 VGM：

```powershell
python vtk_selection_to_vgm.py `
  --occupancy "output\test2_stp_occupancy_target100000000_529x381x496.vtk" `
  --constraint-selection "constraint_selection.vtk" `
  --load-selection "load_selection.vtk" `
  --output "output\test2_stp.vgm" `
  --fine 5 5 5 `
  --spc-dof 123 `
  --force 1000 `
  --direction 0,0,-1
```

说明：

- `--fine 5 5 5` 是默认值，表示每个 VGM 粗单元包含 `5x5x5` 个细体素。脚本会按 5 个 VTK 体素合并成 1 个 VGM 粗体素，不足整除的边界会自动补空。
- `--fine 1 1 1` 表示一个 VTK 体素直接对应一个 VGM 粗体素，最保真，但亿级体素会生成非常大的 VGM，写出和计算都很慢。
- `--spc-dof 123` 表示约束 X/Y/Z 平动自由度。
- `--spc-node-mode cell-corners` 表示约束映射到命中粗单元的 8 个角点；`cell-min-corner` 表示只映射到命中粗单元的最小角点。
- `--force` 是 `FORCE` 卡片中的力大小，`--direction` 是方向向量。

如果约束和载荷来自 UG 导出的两个选中面 STL，可以用：

```powershell
python stl_surface_to_vgm.py `
  --occupancy "output\test2_stp_occupancy_target100000000_529x381x496.vtk" `
  --constraint-stl "..\test2_stp_1.stl" `
  --load-stl "..\test2_stp_2.stl" `
  --output "output\test2_stp_surface_bc_fine5.vgm" `
  --fine 5 5 5 `
  --spc-dof 123 `
  --spc-node-mode cell-corners `
  --force 1000 `
  --direction 0,0,-1 `
  --workers 8
```

## DLUTFEM 计算 VGM

`DLUTFEM.exe` 对中文路径不稳定，建议统一用脚本运行。脚本会把求解器和输入 VGM 临时复制到 `C:\temp\dlutfem_ascii`，计算完成后把 `.h5/.vtk/.db/.log` 回拷到输出目录。

```powershell
python run_dlutfem.py `
  --vgm "output\test2_stp_surface_bc_fine5.vgm" `
  --output-dir "output\dlutfem_test2" `
  --job test2_surface_bc
```

如果第一次运行或求解器目录有更新，可以加 `--refresh-bin`：

```powershell
python run_dlutfem.py --vgm "output\test2_stp_surface_bc_fine5.vgm" --refresh-bin
```

我已经用官方示例和当前生成器写出的小模型验证过：

- 官方示例 `Cube_3x1x1Hexa8Ele_PIML_new(gai)_1119.vgm`：返回码 `0`，导出 `.h5` 和 `.vtk`。
- 当前生成器小模型 `15x5x5, fine=5x5x5`：返回码 `0`，导出 `.h5` 和 `.vtk`。

求解器日志中可能出现 `No element nodes found for the specified type.` 这类信息，目前官方示例也会出现，但不影响静力求解完成和结果导出。

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
