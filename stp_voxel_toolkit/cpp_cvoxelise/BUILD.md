# Build cvoxelise.pyd

这里保存的是 `pyvoxel` 的 C++/pybind11 体素化核心源码。

## 已包含文件

- `Voxel.cpp/.h`：体素化核心。
- `MeshReader.cpp/.h`：STL 读取。
- `MeshWriter.cpp/.h`：网格写出辅助。
- `bind.cpp`：Python 扩展绑定。
- `cvoxelise.vcxproj`：Visual Studio C++ 工程。

## 编译环境

当前工程按本机环境配置：

- Python：`C:\Python314`
- pybind11：`C:\Python314\Lib\site-packages\pybind11\include`
- Eigen：`..\third_party\eigen`
- Visual Studio Build Tools 2019，平台工具集 `v142`

如果你的 Python 路径不同，需要修改 `cvoxelise.vcxproj` 中：

- `AdditionalIncludeDirectories`
- `AdditionalLibraryDirectories`
- `AdditionalDependencies` 中的 `python314.lib`

## 编译命令

本机也可使用 MinGW-w64（`g++` 在 PATH 中）重新生成仓库内的 Python 3.14 运行库。从仓库根目录执行：

```powershell
.\.venv\Scripts\python.exe -m pip install pybind11
& .\stp_voxel_toolkit\cpp_cvoxelise\build_mingw.ps1
```

脚本从 `.venv` 自动获取 Python 和 pybind11 路径，静态链接 MinGW 支持库，并将构建结果复制到 `pyvoxel_runtime/cvoxelise.pyd`。该构建使用串行核心；工具的默认调用本来就使用 `parallel=False`。
新增的 `bounds` 参数需要新版扩展；仅更新 Python 文件而沿用旧版 `.pyd` 不可用。仓库运行库已同步重新编译。

从工具包根目录运行：

```powershell
cd "D:\File\Project\Section\体素化代码\stp_voxel_toolkit\cpp_cvoxelise"
& "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\MSBuild\Current\Bin\MSBuild.exe" cvoxelise.vcxproj /p:Configuration=Release /p:Platform=x64
```

编译成功后会生成：

```text
build\Release\cvoxelise.pyd
```

把它复制到运行时目录：

```powershell
Copy-Item -Force build\Release\cvoxelise.pyd ..\pyvoxel_runtime\cvoxelise.pyd
```

## 运行注意

工具中的调用方式固定为：

```python
grid = voxel.Voxel(nx, ny, nz).VOXELISE(mesh_vertices, ray="xyz", parallel=False)
```

指定网格包围盒时，传入 `bounds=np.array([最小角点, 最大角点])`。包围盒必须完整包含网格；采样点为 `min + (index + 0.5) * (max - min) / count`，与输出 VTK/VGM 的单元间距一致。

当前环境测试发现 `parallel=True` 会崩溃，所以默认禁用并行分支。
