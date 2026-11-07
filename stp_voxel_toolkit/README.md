# STP Voxel VTK Toolkit

这个文件夹是整理后的可用版本，包含 GUI、命令行工具、时间对比脚本、C++ 源码、已编译 C++ 运行库。输入支持 `.stp/.step/.stl`。

## 目录

- `gui.py`：图形界面，选择 `.stp/.step/.stl` 后输出体素 VTK 或直接生成 VGM。
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

源码版要求 **64 位 Python 3.14**。仓库自带的 `cvoxelise.pyd` 是为 Python 3.14 编译的，不能由 Python 3.12/3.13 加载。

首次使用时，直接双击项目根目录的 `setup_env.bat`。它会创建根目录下的 `.venv`、升级 pip、安装依赖并验证 C++ 体素化扩展。

也可以在 PowerShell 中手动执行：

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r stp_voxel_toolkit\requirements.txt
```

安装依赖后双击 `stp_voxel_toolkit\run_gui.bat`。启动脚本会优先使用根目录的 `.venv`，环境不存在或不完整时会自动调用一键安装脚本。

如果使用源码版压缩包，解压后先安装上面的 Python 依赖，再双击 `run_stp_voxel_toolkit.bat`。

如果使用便携版压缩包 `stp_voxel_toolkit_portable.zip`，解压后直接双击 `run_stp_voxel_toolkit.bat`。便携版自带 Python、numpy、gmsh、joblib、体素化运行库和 `solver_bin\DLUTFEM.exe`，不需要在对方电脑安装 Python。

如果给外部用户分发，优先使用安装器 `STPVoxelToolkit_Setup.exe`：双击后会安装到当前用户的 `%LOCALAPPDATA%\STPVoxelToolkit`，并创建桌面和开始菜单快捷方式，不需要管理员权限。

## GUI 使用

```powershell
..\.venv\Scripts\python.exe gui.py
```

GUI 有两个分页：

- `模型生成`：统一选择 STEP/STL 或 occupancy VTK，设置输出目录后点击“开始生成”。默认输出 VGM，也可为 CAD 输入额外选择占据场或表面 VTK。VGM 自动命名为 `<输入文件名>.vgm`。
- `有限元求解`：读取 VGM 并调用 DLUTFEM 计算，结果写入指定目录。

CAD 输入可选择“目标体素数”或“最长轴格数”，两种方式只显示一个数值框。例如目标体素数 `100000000` / `1e8` 对应约一亿总体素。
occupancy VTK 输入沿用原网格，自动隐藏体素化和包围盒设置。
射线、C++ 后端、复用 STL、曲面网格尺寸、选区映射，以及求解器和临时目录等技术选项收在“高级设置”中；默认已启用 C++ 和复用 STL。
网格统计只在结果区显示，模块耗时和运行进度保留在日志中，完成时不再重复弹窗。

### 直接生成 VGM、自定义包围盒和统计

- 在“模型生成”页选择 STEP/STL，保留默认的 VGM 输出即可一次生成；两种 VTK 输出按需勾选。
- 已有 occupancy VTK 时直接将其作为输入。约束和载荷选区对两类输入都分别可选；只有填写选区后，才显示对应的自由度或力参数。
- CAD 输入可勾选“自定义包围盒”，填写世界坐标下的“最小角点 X,Y,Z”和“尺寸 X,Y,Z”，单位与模型相同。模型不移动；包围盒未完全包含模型时立即报错，并列出模型和包围盒范围。尺寸必须为正且所有数值有限。
- 目标体素数和最长轴分辨率按所选包围盒计算。自定义包围盒直接生成 VGM 时，各轴网格数向上对齐至细分的整数倍，再计算间距，因此最终 BOX 保持指定尺寸；实际网格数会在日志里显示。
- 每次生成 VGM 后显示统计，并在同目录保存同名 `*.stats.json`：实体细单元、实体去重节点、细尺度自由度、有效粗单元/节点/自由度，以及含空单元的完整粗网格规模。
- 节点按实体单元的八角点并集去重；自由度按每节点 3 个平移分量统计，是扣除约束及求解器凝聚前的规模。约束粗节点及其位移自由度另列，不能直接从细尺度自由度里相减。

命令行直接输出 VGM（从工具包目录执行，包围盒参数需按实际模型调整）：

```powershell
python stp_to_voxel_vtk.py --input model.stp --resolution 100 `
  --box-origin -10 -10 -10 --box-size 200 200 200 `
  --vgm output/model.vgm --fine-preset 555 --no-occupancy
```

已有 occupancy VTK，不添加约束和载荷：

```powershell
python vtk_selection_to_vgm.py --occupancy output/model_occupancy.vtk --output output/model.vgm
```

VTK 使用世界坐标；VGM 的 `BOX` 定义局部网格尺寸，原点保存在 `$ ORIGIN` 注释和统计 JSON 中。模型在包围盒内的相对位置保持一致。仅含几何的 VGM 可用于后续添加边界条件；静力求解前仍需设置合适的约束和载荷。

“粗单元细分”可选择 `5x5x5`（默认）或 `10x10x10`。只有选择“自定义”时才显示细分 X/Y/Z 输入框。
`10x10x10` 表示每个粗单元含 1000 个细体素，不是整个模型只划分 10 格。
使用同一份 occupancy VTK 时，细体素间距和实体细体素数量保持不变，仅改变粗单元分组。
不足整除的边界补空，BOX 外包尺寸相应增大，不会压缩原有几何。
粗节点及其约束映射会随粗网格改变，因此不同细分的结果不一定完全相同。

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
- `--fine-preset 101010` 或 `--fine 10 10 10`：每个粗单元含 `10x10x10` 个细体素，MATB 和 FORCE 的每个粗单元掩码为 1000 项。
- `--fine-preset 555` 等价于 `--fine 5 5 5`。`--fine-preset` 与 `--fine` 不能同时使用；这两个入口脚本支持相同的细分选项。
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

生成 10×10×10 时，上面任一命令中的 `--fine 5 5 5` 替换为：

```powershell
--fine-preset 101010
```

必须使用支持 10×10×10 的 SGSim/PIML 版本，并随运行时提供 `Data/EN_10` 权重。
生成器不会升级或替换求解器；旧便携包内的 DLUTFEM 不一定支持该细分。
保留约 1 亿细尺度自由度时，应沿用原细网格分辨率，不要因为从 5 改成 10 再加密体素。
`--target-cells` 控制包围盒总体素数，并不是实体节点数或自由度数。

### 细分生成回归测试

从仓库根目录运行：

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s stp_voxel_toolkit/tests -v
```

测试覆盖两种 VGM 写出方式、VTK/STL 选区命令行、GUI 预设、5/10 细分和自定义细分，
并从输出文件还原材料及载荷掩码，检查边界补空、细体素间距和约束节点。

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
