@echo off
setlocal
cd /d "%~dp0"

set "NO_PAUSE="
if /i "%~1"=="--no-pause" set "NO_PAUSE=1"

set "VENV_DIR=%~dp0.venv"
set "VENV_PYTHON=%VENV_DIR%\Scripts\python.exe"
set "REQUIREMENTS=%~dp0stp_voxel_toolkit\requirements.txt"

echo [1/4] Checking 64-bit Python 3.14...
where py >nul 2>nul
if errorlevel 1 goto :python_error

py -3.14 -c "import sys; raise SystemExit(0 if sys.maxsize.bit_length() == 63 else 1)" >nul 2>nul
if errorlevel 1 goto :python_error

echo [2/4] Creating local environment: .venv
if not exist "%VENV_PYTHON%" (
    py -3.14 -m venv "%VENV_DIR%"
    if errorlevel 1 goto :failed
)

"%VENV_PYTHON%" -c "import sys; raise SystemExit(0 if sys.version_info[:2] == (3, 14) and sys.maxsize.bit_length() == 63 else 1)" >nul 2>nul
if errorlevel 1 (
    echo [ERROR] The existing .venv was not created with 64-bit Python 3.14.
    echo Rename or remove .venv, then run this file again.
    goto :failed
)

echo [3/4] Updating pip and installing dependencies...
"%VENV_PYTHON%" -m pip install --upgrade pip
if errorlevel 1 goto :failed
"%VENV_PYTHON%" -m pip install -r "%REQUIREMENTS%"
if errorlevel 1 goto :failed

echo [4/4] Verifying Python packages and the C++ voxelizer...
"%VENV_PYTHON%" -c "import sys; sys.path.insert(0, r'%~dp0stp_voxel_toolkit\pyvoxel_runtime'); import numpy, gmsh, joblib, cvoxelise; print('Environment OK: Python', sys.version.split()[0], '| NumPy', numpy.__version__, '| Gmsh', gmsh.__version__)"
if errorlevel 1 goto :failed

echo.
echo Setup completed successfully.
echo Start the application with stp_voxel_toolkit\run_gui.bat
if not defined NO_PAUSE pause
exit /b 0

:python_error
echo [ERROR] 64-bit Python 3.14 was not found.
echo Download it from https://www.python.org/downloads/windows/
echo During installation, enable the Python Launcher option.

:failed
echo.
echo Environment setup failed.
if not defined NO_PAUSE pause
exit /b 1
