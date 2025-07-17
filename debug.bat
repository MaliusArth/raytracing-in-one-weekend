@echo off

echo %~0 %*

:: setlocal enabledelayedexpansion
if "%1" == "" (
    echo missing source file
    goto :eof
)

set "project_path=%~dp0"
set "exe_name=%~n1"
set "source_path=%~2"
:: set "source_name=%~n2"
if "%source_path%"=="" set "source_name=%exe_name%"


:: tasklist /fi "imagename eq raddbg.exe"

for /f %%i in ('qprocess "remedybg.exe" ^2^>nul ^| find /i /c "remedybg.exe"') do set "runningProcesses=%%i"
:: echo running debugger processes: %runningProcesses%
:: echo project path: %project_path%
if "%runningProcesses%"=="0" (
    start /max remedybg "bin\%exe_name%.exe"
) else (
    remedybg bring-debugger-to-foreground
)

:: for %%A in (*.rdbg) do set session_file=%%~fsA
:: if exist "%session_file%" start remedybg open-session "%session_file%"

:: set /p input=press ENTER to continue...
:: start remedybg "bin\%exe_name%.exe"
:: start remedybg open-file "bin\%exe_name%.exe"
:: timeout /t 1
:: start /wait remedybg bring-debugger-to-foreground
remedybg stop-debugging
remedybg open-file "%source_path%"
remedybg add-breakpoint-at-function "main.main"
call remedybg start-debugging
:: remedybg remove-breakpoint-at-function "main.main"



:: TODO: else activate program window & exec ipc call (maybe the ipc call also activates it?)

exit /b
:: set




:: @if (@a==@a) @end /*
::     @echo off
::     cscript //E:JScript //nologo "%~f0" %*

::     exit /b %errorlevel%
:: */

:: // --- JScript code below this line ----

:: var WshShell = WScript.CreateObject("WScript.Shell");
:: // var ARGS = WScript.Arguments;

:: WshShell.AppActivate("The RAD Debugger");
:: // WshShell.SendKeys("{F5}");
:: // WScript.Echo(ARGS.Item(0) + " activated");
:: WScript.Quit(0);