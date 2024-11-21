@echo off
@REM TODO: update to properly support package folders

echo %~0 %*

@REM setlocal enabledelayedexpansion
if "%1" == "" (
    echo "you forgot to pass a source file"
    goto :eof
)

set "project_path=%~dp0"
set "exe_name=%~n1"
set "source_path=%~2"
@REM set "source_name=%~n2"
if "%source_path%"=="" set "source_name=%exe_name%"


@REM tasklist /fi "imagename eq raddbg.exe"

for /f %%i in ('qprocess "remedybg.exe" ^2^>nul ^| find /i /c "remedybg.exe"') do set "runningProcesses=%%i"
@REM echo running debugger processes: %runningProcesses%
@REM echo project path: %project_path%
if "%runningProcesses%"=="0" (
    start /max remedybg "bin\%exe_name%.exe"
) else (
    remedybg bring-debugger-to-foreground
)

@REM for %%A in (*.rdbg) do set session_file=%%~fsA
@REM if exist "%session_file%" start remedybg open-session "%session_file%"

@REM set /p input=press ENTER to continue...
@REM start remedybg "bin\%exe_name%.exe"
@REM start remedybg open-file "bin\%exe_name%.exe"
@REM timeout /t 1
@REM start /wait remedybg bring-debugger-to-foreground
remedybg stop-debugging
remedybg open-file "%source_path%"
remedybg add-breakpoint-at-function "main.main"
remedybg start-debugging
@REM remedybg remove-breakpoint-at-function "main.main"



@REM TODO: else activate program window & exec ipc call (maybe the ipc call also activates it?)

exit /b
@REM set




@REM @if (@a==@a) @end /*
@REM     @echo off
@REM     cscript //E:JScript //nologo "%~f0" %*

@REM     exit /b %errorlevel%
@REM */

@REM // --- JScript code below this line ----

@REM var WshShell = WScript.CreateObject("WScript.Shell");
@REM // var ARGS = WScript.Arguments;

@REM WshShell.AppActivate("The RAD Debugger");
@REM // WshShell.SendKeys("{F5}");
@REM // WScript.Echo(ARGS.Item(0) + " activated");
@REM WScript.Quit(0);