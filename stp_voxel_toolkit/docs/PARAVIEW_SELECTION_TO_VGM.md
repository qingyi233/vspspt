# ParaView 选区生成 VGM 使用说明

本文说明如何用 ParaView 在体素表面 VTK 上选择约束区域和载荷区域，然后用脚本生成自定义 VGM 文件。

## 一、整体流程

```text
STL / STP
  -> 体素化，生成 occupancy VTK 和 surface VTK
  -> ParaView 打开 surface VTK
  -> 选择约束区域，保存 constraint_selection.vtk
  -> 选择载荷区域，保存 load_selection.vtk
  -> Python 脚本读取 occupancy VTK + 两个选区 VTK
  -> 输出 VGM
```

注意：不能只靠两个选区 VTK 生成完整 VGM。完整几何、材料分布、体素占据关系来自 `*_occupancy_*.vtk`；两个选区 VTK 只负责生成约束和载荷。

## 二、重新生成带编号的 surface VTK

先进入工具目录：

```powershell
cd "D:\File\Project\Section\体素化代码\stp_voxel_toolkit"
```

如果输入是 STL，例如：

```powershell
python stp_to_voxel_vtk.py --input "..\test2_stp.stl" --target-cells 100000000 --surface
```

如果输入是 STP，并且已经有可复用 STL：

```powershell
python stp_to_voxel_vtk.py --input "..\test1.stp" --target-cells 100000000 --surface --reuse-stl
```

输出目录一般在：

```text
D:\File\Project\Section\体素化代码\stp_voxel_toolkit\output
```

会得到两类文件：

```text
xxx_occupancy_xxx.vtk
xxx_surface_xxx.vtk
```

其中：

- `*_occupancy_*.vtk`：完整体素占据场，用于生成 `BOX / COARSE / MATB`。
- `*_surface_*.vtk`：给 ParaView 显示和选择用。

新版 `*_surface_*.vtk` 内部带有这些 cell data：

```text
voxel_i
voxel_j
voxel_k
voxel_eid
face_id
```

这些字段用于把 ParaView 选中的表面区域映射回 VGM 体素编号。旧版本生成的 surface VTK 没有这些字段，需要重新生成。

## 三、在 ParaView 中选择约束区域

1. 打开 ParaView。
2. `File -> Open`，选择生成的 `*_surface_*.vtk`。
3. 点击 `Apply`。
4. 用 ParaView 的选择工具选择要施加约束的区域。
5. 选择完成后，使用：

```text
Filters -> Data Analysis -> Extract Selection
```

6. 点击 `Apply`。
7. 选中左侧 Pipeline 中的 `ExtractSelection` 结果。
8. 使用：

```text
File -> Save Data
```

9. 保存为 legacy VTK 文件，例如：

```text
constraint_selection.vtk
```

建议保存到工具目录，或者保存到你方便写路径的位置。

## 四、在 ParaView 中选择载荷区域

重复第三节的操作，只是这次选择要施加载荷的区域。

保存为：

```text
load_selection.vtk
```

最终你应该有三个输入文件：

```text
xxx_occupancy_xxx.vtk
constraint_selection.vtk
load_selection.vtk
```

## 五、生成 VGM

### 方法 A：在 GUI 中生成

运行 GUI：

```powershell
python gui.py
```

在界面中的 `ParaView 选区生成 VGM` 区域填写：

- `占据场 VTK`：选择原始 `*_occupancy_*.vtk`。
- `约束选区 VTK`：选择 ParaView 保存的 `constraint_selection.vtk`。
- `载荷选区 VTK`：选择 ParaView 保存的 `load_selection.vtk`。
- `输出 VGM`：选择要保存的 `.vgm` 文件。
- `细分`：第一轮建议填 `1 1 1`。
- `约束自由度`：常用 `123`。
- `力大小`：例如 `1000`。
- `方向`：例如 `0,0,-1`。

点击 `生成 VGM`。

如果你先用 GUI 完成了体素化，生成的占据场 VTK 会自动填入 `占据场 VTK` 输入框。你只需要在 ParaView 中保存两个选区，然后回到 GUI 选择这两个选区文件。

GUI 参数含义：

- `细分 X/Y/Z`：VGM 中一个粗体素内部再划分的小体素数量。`1 1 1` 表示一个 VTK 体素就是一个 VGM 粗体素；`5 5 5` 表示每个 VGM 粗体素内部有 125 个细体素。
- `约束自由度`：写入 `SPC1` 的自由度。`123` 表示固定 X/Y/Z 三个平动自由度。
- `SPC粗节点映射`：现有求解器的 `SPC1` 只支持粗网格节点，不支持细体素内部节点。`命中粗单元8角点` 表示凡是被约束曲面命中的粗单元，其 8 个粗节点都写入 `SPC1`；`命中粗单元最小角点` 只写入每个命中粗单元的最小角点，约束更少。
- `力大小`：写入 `FORCE` 卡片的力幅值。
- `方向`：写入 `FORCE` 卡片的方向向量，例如 `0,0,-1` 表示 Z 轴负方向。

### 方法 B：命令行生成

进入工具目录：

```powershell
cd "D:\File\Project\Section\体素化代码\stp_voxel_toolkit"
```

运行：

```powershell
python vtk_selection_to_vgm.py `
  --occupancy "output\xxx_occupancy_xxx.vtk" `
  --constraint-selection "constraint_selection.vtk" `
  --load-selection "load_selection.vtk" `
  --output "output\model.vgm" `
  --fine 1 1 1 `
  --spc-dof 123 `
  --force 1000 `
  --direction 0,0,-1
```

把 `output\xxx_occupancy_xxx.vtk` 换成你真实生成的 occupancy 文件名。

## 六、参数说明

### `--occupancy`

原始体素占据文件，例如：

```text
output\test2_stp_occupancy_target100000000_529x381x496.vtk
```

这个文件用于生成 VGM 的几何和材料体素分布。

### `--constraint-selection`

ParaView 中保存的约束选区 VTK。

该选区会转换为：

```text
SPCADD
SPC1
```

### `--load-selection`

ParaView 中保存的载荷选区 VTK。

该选区会转换为：

```text
LOAD
FORCE
```

### `--fine`

VGM 中每个粗体素内部的细分数量。

推荐第一轮测试用：

```powershell
--fine 1 1 1
```

含义是：一个 VTK 体素直接对应一个 VGM 粗体素，最不容易出错。

如果你必须生成示例文件那种 `5x5x5` 细分格式，可以用：

```powershell
--fine 5 5 5
```

这时脚本会把 `5x5x5` 个 VTK 体素合成一个 VGM 粗体素。不足整除的边界会自动补空。

### `--spc-dof`

约束自由度。

常用：

```text
123
```

表示约束 X/Y/Z 三个平动自由度。

### `--force`

载荷大小，也就是 VGM `FORCE` 卡片中的 `F`。

例如：

```powershell
--force 1000
```

### `--direction`

载荷方向向量。

例如 Z 轴负方向：

```powershell
--direction 0,0,-1
```

X 轴正方向：

```powershell
--direction 1,0,0
```

## 七、推荐测试命令

先用较小体素数测试完整流程，例如：

```powershell
python stp_to_voxel_vtk.py --input "..\test2_stp.stl" --target-cells 1000000 --surface
```

在 ParaView 中打开 `*_surface_*.vtk`，分别保存：

```text
constraint_selection.vtk
load_selection.vtk
```

然后生成 VGM：

```powershell
python vtk_selection_to_vgm.py `
  --occupancy "output\test2_stp_occupancy_target1000000_xxx.vtk" `
  --constraint-selection "constraint_selection.vtk" `
  --load-selection "load_selection.vtk" `
  --output "output\test2_stp_test.vgm" `
  --fine 1 1 1 `
  --spc-dof 123 `
  --force 1000 `
  --direction 0,0,-1
```

其中 `xxx` 部分以实际输出文件名为准。

## 八、常见问题

### 1. 报错说选区没有 `voxel_i / voxel_j / voxel_k`

原因通常是你用的是旧版生成的 `*_surface_*.vtk`。

解决办法：重新运行：

```powershell
python stp_to_voxel_vtk.py --input "你的模型.stl" --target-cells 1000000 --surface
```

然后在 ParaView 中重新选择并保存选区。

### 2. ParaView 中看不到 occupancy VTK

`*_occupancy_*.vtk` 是结构化体素场，不一定直接显示出实体。

ParaView 中可以对它做 Threshold，选择：

```text
occupancy = 1
```

不过选区建议使用 `*_surface_*.vtk`，它更直观。

### 3. VGM 文件太大

体素数越多，`COARSE` 和 `MATB` 越多，VGM 会非常大。

测试流程时建议先用：

```powershell
--target-cells 1000000
```

流程确认没问题后再提高到：

```powershell
--target-cells 100000000
```

### 4. `--fine 5 5 5` 后粗体素数量变化

这是正常的。

例如原始 VTK 是 `529 x 381 x 496`，用 `--fine 5 5 5` 后，VGM 粗体素数量会变成：

```text
ceil(529/5) x ceil(381/5) x ceil(496/5)
```

边界不足 5 的部分会补空。
