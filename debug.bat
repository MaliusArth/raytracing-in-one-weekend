@echo off
@REM setlocal enabledelayedexpansion
if "%1" == "" (
    echo "you forgot to pass a source file"
    goto :eof
)

set "file_name=%~n1"
set "project_path=%~dp0"

@REM tasklist /fi "imagename eq raddbg.exe"

for /f %%i in ('qprocess "remedybg.exe" ^| find /i /c "remedybg.exe"') do set "runningProcesses=%%i"
echo running debugger processes: %runningProcesses%
echo project path: %project_path%
if "%runningProcesses%"=="0" start /max remedybg.exe "%project_path%bin\%file_name%.exe" else remedybg.exe bring-debugger-to-foreground

@REM for %%A in (*.rdbg) do set session_file=%%~fsA
@REM if exist "%session_file%" start remedybg.exe open-session "%session_file%"

@REM set /p input=press ENTER to continue...
@REM start remedybg.exe "bin\%file_name%.exe"
@REM start remedybg.exe open-file "bin\%file_name%.exe"
@REM timeout /t 1
@REM start /wait remedybg.exe bring-debugger-to-foreground
remedybg.exe stop-debugging
remedybg.exe open-file "%project_path%src\%file_name%.odin"
remedybg.exe add-breakpoint-at-function "main.main"
remedybg.exe start-debugging
@REM remedybg.exe remove-breakpoint-at-function "main.main"



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