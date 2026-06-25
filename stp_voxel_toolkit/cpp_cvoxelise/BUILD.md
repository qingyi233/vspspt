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

当前环境测试发现 `parallel=True` 会崩溃，所以默认禁用并行分支。
