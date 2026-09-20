@echo off
setlocal
cd /d "%~dp0"

if exist "%~dp0python\pythonw.exe" (
    start "" "%~dp0python\pythonw.exe" "%~dp0gui.py"
    exit /b 0
)

call "%~dp0run_gui.bat"
