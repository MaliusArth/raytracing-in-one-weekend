@echo off

:: batch comments ref
:: https://stackoverflow.com/questions/12407800/which-comment-style-should-i-use-in-batch-files
:: https://stackoverflow.com/questions/64186507/question-about-comments-in-batch-bat-files-and-speed/68436626#68436626

echo %~0 %*

setlocal EnableDelayedExpansion

:: usage

:: build path/to/file -file -test
:: build path/to/package -test

:: defaults

:: set "-test=0"
set "source=%~1"
:: if not exist %source% echo source '%source%' does not exist! && goto:eof
shift

:parseArgs
call:getArgFlag "-file" "-file" "%~1" && shift && goto :parseArgs
call:getArgFlag "-test" "-test" "%~1" && shift && goto :parseArgs
:: call:getArgWithValue "-file" "-file" "%~1" "%~2" && shift && shift && goto :parseArgs

:: echo -file: %-file%
:: echo -test: %-test%

if "%-file%" == "1" (
    for %%A in ("%source%\.") do set "out_name=%%~nA"
    set "source=%source% -file"
) else (
    set "out_name=%source%"
)

set "command=build"
if "%-test%" == "1" set "command=test"

:: set "project_path=%~dp0"
:: for %%A in ("%~p0\.") do set "project_name=%%~nxA"

set "out_dir=bin"

:: -vet
::   -vet-unused
::   -vet-unused-variables
::   -vet-unused-imports
::   -vet-shadowing
::   -vet-using-stmt

:: TODO: use in git push hook to prevent 'using' getting committed to the repo
:: -vet-using-stmt & -vet-using-param
::   'using' is considered bad practice outside of immediate refactoring.

:: -vet-style
:: -vet-semicolon
:: -vet-cast
:: -vet-tabs

:: -strict-style
::   like -vet-style -vet-semicolon
::   + Errs when the attached-brace style in not adhered to (also known as 1TBS).
::   + Errs when 'case' labels are not in the same column as the associated 'switch' token.
:: -vet-style -vet-semicolon
set STYLE_PARAMS=-vet -strict-style -vet-using-param -vet-cast -vet-tabs

:: TODO: -sanitize:address

:: -vet
:: odin build src/%1.odin -debug -show-timings -out=bin/%1.exe -extra-linker-flags="/libpath:lib"

:: copy "lib\Debug\raylib_static.pdb" "bin\" >nul 2>&1
:: odin build src/%1.odin -debug -show-timings -out=bin/%1.exe -extra-linker-flags="/libpath:lib\Debug /NODEFAULTLIB:raylib_static" -subsystem:windows -keep-temp-files

:: SET "raylib_dir=lib\raylib"
:: copy "%raylib_dir%\raylib.pdb" "bin\" >nul 2>&1
:: copy "%raylib_dir%\raylib.dll" "bin\" >nul 2>&1
::  -keep-temp-files
:: -define:RAYLIB_SHARED=true to link against a dll version
:: -define:RAYLIB_SHARED=true
:: odin build src/%1.odin -file -out=bin/%1.exe -o:none -debug -show-timings -extra-linker-flags="/libpath:%raylib_dir%" -subsystem:console %STYLE_PARAMS%
:: odin strip-semicolon src/%file_name%.odin -file

mkdir %out_dir% 2>nul

:: -o:none
:: -o:minimal (default)
:: -o:size
:: -o:speed
:: -o:aggressive
odin %command% %source% -out=%out_dir%/%out_name%.exe -o:speed -debug -show-timings -subsystem:console %STYLE_PARAMS%

:: start rundll32.exe cmdext.dll,MessageBeepStub
start rundll32 user32.dll,MessageBeep
:: call powershell "[console]::beep(500,300)"

:: || goto :error
:: echo -------------------------

:: goto :end
:: :error
:: if %errorlevel% neq 0 exit /b %errorlevel%

:: :end



goto:eof
:: ===== END ===== ::

:: "-foo" "FOO" "%~1" "%~2"

:: =====================================================================
:: Sets a variable from a cli argument with value
:: 1 variable name
:: 2 parameter name
:: 3 argument name
:: 4 argument value
:: 5 default value
::   if "%~4"=="" (
::     REM unset the variable if value is not provided
::     set "%~1="
::     exit /B 0
::   )
:getArgWithValue
if "%~3"=="%~2" (
  set "%~1=%~4"
  exit /B 0
)
exit /B 1

:: =====================================================================
:: Sets a variable from a cli flag argument
:: 1 cli argument name
:: 2 variable name
:: 3 current Argument Name
:getArgFlag
if "%~3"=="%~2" (
  set "%~1=1"
  exit /B 0
)
exit /B 1