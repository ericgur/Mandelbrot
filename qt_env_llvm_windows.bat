@echo off

set QT_ROOT=C:\Qt
set QT_VER=
set QT_TOOLCHAIN=
setlocal EnableDelayedExpansion

set "latestVersionNum=000000000"
set "latestVersion="

:: Loop through all subdirectories
for /d %%D in ("%QT_ROOT%\*") do (
    set "dirName=%%~nxD"

    :: Split into major.minor.patch using tokens
    for /f "tokens=1-3 delims=." %%a in ("!dirName!") do (
        :: Ensure all tokens are numbers and third token exists
        echo %%a %%b %%c | findstr /r "^[0-9][0-9]* [0-9][0-9]* [0-9][0-9]*$" >nul
        if !errorlevel! == 0 (
            :: Pad each number to 3 digits
            set "maj=00%%a"
            set "min=00%%b"
            set "pat=00%%c"
            set "verNum=!maj:~-3!!min:~-3!!pat:~-3!"

            :: Compare with latest
            if "!verNum!" GTR "!latestVersionNum!" (
                set "latestVersionNum=!verNum!"
                set "QT_VER=!dirName!"
            )
        )
    )
)

:: Find the latest llvm-mingw toolchain installed under %QT_ROOT%\Tools
set "latestToolchainNum=00000000"

for /d %%D in ("%QT_ROOT%\Tools\llvm-mingw*_64") do (
    set "dirName=%%~nxD"

    :: Strip the "llvm-mingw" prefix, then drop the "_64" suffix
    set "tcVer=!dirName:~10!"
    for /f "tokens=1 delims=_" %%a in ("!tcVer!") do set "tcVer=%%a"

    :: Ensure the version part is numeric
    echo !tcVer!| findstr /r "^[0-9][0-9]*$" >nul
    if !errorlevel! == 0 (
        :: Pad to 8 digits so the string compare below orders numerically
        set "tcNum=0000000!tcVer!"

        :: Compare with latest
        if "!tcNum:~-8!" GTR "!latestToolchainNum!" (
            set "latestToolchainNum=!tcNum:~-8!"
            set "QT_TOOLCHAIN=!dirName!"
        )
    )
)

:: cleanup
set maj=
set min=
set pat=
set verNum=
set latestVersionNum=
set dirName=
set tcVer=
set tcNum=
set latestToolchainNum=

if not defined QT_VER (
    :: failed to find a version
    echo No valid Qt versions found in %QT_ROOT%
    pause
    exit /b 1
)

if not defined QT_TOOLCHAIN (
    :: failed to find a toolchain
    echo No llvm-mingw toolchain found in %QT_ROOT%\Tools
    pause
    exit /b 1
)

echo Latest Qt version found: %QT_VER%
echo Latest llvm-mingw toolchain found: %QT_TOOLCHAIN%
echo Setting up environment for Qt %QT_VER% using toolchain %QT_TOOLCHAIN%...
title Build env for %QT_VER% using toolchain %QT_TOOLCHAIN%
prompt $CQT$F $P$G

:: set the environment variables
cmd /K set PATH=%QT_ROOT%\Tools\Ninja;%QT_ROOT%\Tools\CMake_64;%QT_ROOT%\%QT_VER%\llvm-mingw_64\bin;%QT_ROOT%\Tools\%QT_TOOLCHAIN%\bin;%PATH%
