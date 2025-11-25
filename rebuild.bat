@echo off
call clean.bat %*
IF %ERRORLEVEL% NEQ 0 (Echo An error was found &Exit /b 1)
call build.bat %*
