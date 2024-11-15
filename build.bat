@echo off

echo %~0 %*

setlocal EnableDelayedExpansion

:: usage

:: build path/to/file -file -test
:: build path/to/package -test

:: defaults

@REM set "-test=0"
set "source=%~1"
@REM if not exist %source% echo source '%source%' does not exist! && goto:eof
shift

:parseArgs
call:getArgFlag "-file" "-file" "%~1" && shift && goto :parseArgs
call:getArgFlag "-test" "-test" "%~1" && shift && goto :parseArgs
@REM call:getArgWithValue "-file" "-file" "%~1" "%~2" && shift && shift && goto :parseArgs

@REM echo -file: %-file%
@REM echo -test: %-test%

if "%-file%" == "1" (
    for %%A in ("%source%\.") do set "out_name=%%~nA"
    set "source=%source% -file"
) else (
    set "out_name=%source%"
)

set "command=build"
if "%-test%" == "1" set "command=test"

@REM set "project_path=%~dp0"
@REM for %%A in ("%~p0\.") do set "project_name=%%~nxA"

set "out_dir=bin"

: -vet
:   -vet-unused
:   -vet-unused-variables
:   -vet-unused-imports
:   -vet-shadowing
:   -vet-using-stmt

: TODO: use in git push hook to prevent 'using' getting commited to the repo
: -vet-using-stmt & -vet-using-param
:   'using' is considered bad practice outside of immediate refactoring.

: -vet-style
: -vet-semicolon
: -vet-cast
: -vet-tabs

: -strict-style
:   like -vet-style -vet-semicolon
:   + Errs when the attached-brace style in not adhered to (also known as 1TBS).
:   + Errs when 'case' labels are not in the same column as the associated 'switch' token.
: -vet-style -vet-semicolon
set STYLE_PARAMS=-vet -strict-style -vet-using-param -vet-cast -vet-tabs

: TODO: -sanitize:address

: -vet
@REM odin build src/%1.odin -debug -show-timings -out=bin/%1.exe -extra-linker-flags="/libpath:lib"

@REM copy "lib\Debug\raylib_static.pdb" "bin\" >nul 2>&1
@REM odin build src/%1.odin -debug -show-timings -out=bin/%1.exe -extra-linker-flags="/libpath:lib\Debug /NODEFAULTLIB:raylib_static" -subsystem:windows -keep-temp-files

@REM SET "raylib_dir=lib\raylib"
@REM copy "%raylib_dir%\raylib.pdb" "bin\" >nul 2>&1
@REM copy "%raylib_dir%\raylib.dll" "bin\" >nul 2>&1
@REM  -keep-temp-files
: -define:RAYLIB_SHARED=true to link against a dll version
@REM -define:RAYLIB_SHARED=true
@REM odin build src/%1.odin -file -out=bin/%1.exe -o:none -debug -show-timings -extra-linker-flags="/libpath:%raylib_dir%" -subsystem:console %STYLE_PARAMS%
@REM odin strip-semicolon src/%file_name%.odin -file

mkdir %out_dir% 2>nul

odin %command% %source% -out=%out_dir%/%out_name%.exe -o:speed -debug -show-timings -subsystem:console %STYLE_PARAMS%

@REM start rundll32.exe cmdext.dll,MessageBeepStub
start rundll32 user32.dll,MessageBeep
@REM call powershell "[console]::beep(500,300)"

@REM || goto :error
@REM echo -------------------------

@REM goto :end
@REM :error
@REM if %errorlevel% neq 0 exit /b %errorlevel%

@REM :end



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
:getArgWithValue
if "%~3"=="%~2" (
@REM   if "%~4"=="" (
@REM     REM unset the variable if value is not provided
@REM     set "%~1="
@REM     exit /B 0
@REM   )
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