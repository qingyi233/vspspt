param(
    [string]$Python = "$PSScriptRoot\..\..\.venv\Scripts\python.exe",
    [string]$Compiler = "g++"
)

$ErrorActionPreference = 'Stop'
$taskPythonRoot = & $Python -c 'import sys; print(sys.base_prefix)'
if ($LASTEXITCODE -ne 0) { throw 'Unable to locate Python' }
$taskPybindInclude = & $Python -c 'import pybind11; print(pybind11.get_include())'
if ($LASTEXITCODE -ne 0) { throw 'Install build dependency: python -m pip install pybind11' }
$taskPythonLib = & $Python -c 'import sys; print("python%d%d.lib" % sys.version_info[:2])'
$taskBuildDir = Join-Path $PSScriptRoot 'build'
New-Item -ItemType Directory -Force -Path $taskBuildDir | Out-Null
$taskOutput = Join-Path $taskBuildDir 'cvoxelise.pyd'
$taskSources = @('bind.cpp', 'Voxel.cpp', 'MeshReader.cpp', 'MeshWriter.cpp', 'base64.cpp') |
    ForEach-Object { Join-Path $PSScriptRoot $_ }
# Statically link MinGW support libraries so the runtime needs no new DLLs.
& $Compiler -shared -std=c++17 -O2 -static-libgcc -static-libstdc++ `
    -I "$taskPythonRoot\include" -I $taskPybindInclude -I "$PSScriptRoot\..\third_party\eigen" `
    @taskSources "$taskPythonRoot\libs\$taskPythonLib" -o $taskOutput `
    '-Wl,-Bstatic' -lwinpthread '-Wl,-Bdynamic'
if ($LASTEXITCODE -ne 0) { throw "C++ build failed: $LASTEXITCODE" }
Copy-Item -LiteralPath $taskOutput -Destination "$PSScriptRoot\..\pyvoxel_runtime\cvoxelise.pyd" -Force
Write-Output 'Built and installed cvoxelise.pyd'
