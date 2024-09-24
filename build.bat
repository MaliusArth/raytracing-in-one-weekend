@echo off
if "%1" == "" (
    echo "you forgot to pass a source file"
    goto :eof
)

set "file_name=%~n1"

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
set STYLE_PARAMS=-vet -vet-style -vet-semicolon -vet-using-param -vet-cast -vet-tabs

: TODO: -sanitize:address

: -vet
@REM odin build src/%1.odin -debug -show-timings -out=bin/%1.exe -extra-linker-flags="/libpath:lib"

@REM copy "lib\Debug\raylib_static.pdb" "bin\" > NUL 2>&1
@REM odin build src/%1.odin -debug -show-timings -out=bin/%1.exe -extra-linker-flags="/libpath:lib\Debug /NODEFAULTLIB:raylib_static" -subsystem:windows -keep-temp-files

@REM SET "raylib_dir=lib\raylib"
@REM copy "%raylib_dir%\raylib.pdb" "bin\" > NUL 2>&1
@REM copy "%raylib_dir%\raylib.dll" "bin\" > NUL 2>&1
@REM  -keep-temp-files
: -define:RAYLIB_SHARED=true to link against a dll version
@REM -define:RAYLIB_SHARED=true
@REM odin build src/%1.odin -file -out=bin/%1.exe -o:none -debug -show-timings -extra-linker-flags="/libpath:%raylib_dir%" -subsystem:console %STYLE_PARAMS%
@REM odin strip-semicolon src/%file_name%.odin -file
odin build src/%file_name%.odin -file -out=bin/%file_name%.exe -o:none -debug -show-timings -subsystem:console %STYLE_PARAMS%

@REM || goto :error
@REM echo -------------------------

@REM goto :end
@REM :error
@REM if %errorlevel% neq 0 exit /b %errorlevel%

@REM :end
