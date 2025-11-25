@echo off
call build.bat %*
IF %ERRORLEVEL% NEQ 0 (Echo An error was found &Exit /b 1)
if "%~1"=="" (
    cmake --install build
) else (
    cmake --install build --config %1
)
