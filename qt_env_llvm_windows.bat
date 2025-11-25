@echo off

set QT_ROOT=C:\Qt
set QT_VER=
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

:: cleanup
set maj=
set min=
set pat=
set verNum=
set latestVersionNum=

if defined QT_VER (
    :: got a version
    echo Latest Qt version found: %QT_VER%
    set QT_TOOLCHAIN=llvm-mingw1706_64
    echo Setting up environment for Qt %QT_VER% using toolchain %QT_TOOLCHAIN%...
    title Build env for %QT_VER% using toolchain %QT_TOOLCHAIN%
    prompt $CQT$F $P$G
) else (
    :: failed to find a version
    echo No valid Qt versions found in %QT_ROOT%
    pause
    exit /b 1
)

:: set the environment variables
cmd /K set PATH=%QT_ROOT%\Tools\Ninja;%QT_ROOT%\Tools\CMake_64;%QT_ROOT%\%QT_VER%\llvm-mingw_64\bin;%QT_ROOT%\Tools\%QT_TOOLCHAIN%\bin;%PATH%
