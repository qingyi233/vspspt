@echo off
setlocal
cd /d "%~dp0"

set "LOCAL_PYTHON=%~dp0..\.venv\Scripts\python.exe"
set "SETUP_SCRIPT=%~dp0..\setup_env.bat"

if not exist "%LOCAL_PYTHON%" (
    echo Local environment not found. Starting one-click setup...
    call "%SETUP_SCRIPT%" --no-pause
    if errorlevel 1 exit /b 1
)

"%LOCAL_PYTHON%" -c "import gui" >nul 2>nul
if errorlevel 1 (
    echo Local environment is incomplete. Repairing it...
    call "%SETUP_SCRIPT%" --no-pause
    if errorlevel 1 exit /b 1
)

"%LOCAL_PYTHON%" gui.py
pause
