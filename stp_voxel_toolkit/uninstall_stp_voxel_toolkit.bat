@echo off
setlocal
set "APP_PARENT=%LOCALAPPDATA%\STPVoxelToolkit"
set "DESKTOP_LINK=%USERPROFILE%\Desktop\STP Voxel Toolkit.lnk"
set "START_DIR=%APPDATA%\Microsoft\Windows\Start Menu\Programs\STP Voxel Toolkit"

echo Removing STP Voxel Toolkit...
if exist "%DESKTOP_LINK%" del /f /q "%DESKTOP_LINK%"
if exist "%START_DIR%" rmdir /s /q "%START_DIR%"
if exist "%APP_PARENT%" rmdir /s /q "%APP_PARENT%"
echo Done.
pause
