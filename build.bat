@echo off
if not exist build (
    call configure.bat
)

if "%~1"=="" (
    cmake --build build
) else (
    cmake --build build --config %1
)
