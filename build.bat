:: https://stackoverflow.com/questions/4094699/how-does-the-windows-command-interpreter-cmd-exe-parse-scripts/4095133#4095133
@echo off
setlocal DisableDelayedExpansion
call :setESC

:: batch comments ref
:: https://stackoverflow.com/questions/12407800/which-comment-style-should-i-use-in-batch-files
:: https://stackoverflow.com/questions/64186507/question-about-comments-in-batch-bat-files-and-speed/68436626#68436626

echo.%~0 %* >&2

:: usage

:: build path/to/file -file -test
:: build path/to/package -test

:: positional args
set "-build_setting=%~1" & shift
set "-source=%~1"  & shift

if "%-source%" == ""     echo missing source package/file        & goto :eof
if not exist "%-source%" echo source '%-source%' does not exist! & goto :eof

:parse_args
:: variable cli_match cur_arg
call :get_arg_flag "-file" "-file" "%~1" && shift && goto :parse_args
call :get_arg_flag "-dry"  "-dry"  "%~1" && shift && goto :parse_args
call :get_arg_flag "-run"  "-run"  "%~1" && shift && goto :parse_args
call :get_arg_flag "-test"  "-test"  "%~1" && shift && goto :parse_args
@REM call :get_arg_with_value "-build_setting" "-build_setting" "%~1" "%~2" "build" && shift && shift && goto :parse_args

@REM echo -file: %-file%
@REM echo -dry: %-dry%
@REM echo -run: %-run%
@REM echo -test: %-test%

:: requires DisableDelayedExpansion
call :define_macro_execute

setlocal EnableDelayedExpansion
:: ============================================================================
:: main script starts here

if "%-file%"=="1" (
    for %%A in ("%-source%\.") do set "out_name=%%~nA"
    set "-source=%-source% -file"
) else (
    set "out_name=%-source%"
)
if "%-test%"=="1" set "out_name=%out_name%_test"

:: -vet
::   -vet-unused
::   -vet-unused-variables
::   -vet-unused-imports
::   -vet-shadowing
::   -vet-using-stmt

:: TODO: use in git push hook to prevent 'using' getting committed to the repo
:: -vet-using-stmt & -vet-using-param
::   'using' is considered bad practice outside of immediate refactoring.
:: odin strip-semicolon src/%file_name%.odin -file

:: -vet-style
:: -vet-semicolon
:: -vet-cast
:: -vet-tabs

:: -strict-style
::   like -vet-style -vet-semicolon
::   + Errs when the attached-brace style in not adhered to (also known as 1TBS).
::   + Errs when 'case' labels are not in the same column as the associated 'switch' token.
:: -vet-style -vet-semicolon
set "style_settings=-vet -strict-style -vet-using-param -vet-cast -vet-tabs"

:: -o:none
:: -o:minimal (default)
:: -o:size
:: -o:speed
:: -o:aggressive

:: -sanitize:address

if "%-build_setting%"=="debug"   set "compiler_settings=-o:minimal -debug"
if "%-build_setting%"=="release" set "compiler_settings=-o:speed   -debug"
if "%-build_setting%"=="fast"    set "compiler_settings=-o:aggressive -disable-assert"
:: if "%-build_setting%"=="test"    set "compiler_settings=-o:aggressive -disable-assert"
if "%-build_setting%"=="asm"     set "compiler_settings=-o:speed -build-mode:asm -keep-temp-files"

set "-command=build"
if "%-test%"=="1" set "-command=test"
if "%-run%"=="1"  set "-command=run"

set "out_dir=bin"
mkdir %out_dir% 2>nul
%execute% odin %-command% %-source% -out=%out_dir%/%out_name%.exe %compiler_settings% -subsystem:console %style_settings%
if %errorlevel%==0 (
  start /b cmd /c call "misc\spplayer.bat" "C:\Windows\Media\tada.wav"
) else (
  start /b cmd /c call "misc\spplayer.bat" "C:\Windows\Media\chord.wav"
)
:: start rundll32 user32.dll,MessageBeep MB_OK
:: start rundll32 user32.dll,MessageBeep MB_ERROR
:: BEL character Alt+007 in cmd line
:: echo 
:: echo lol | clip



exit /b %errorlevel%

:: end
:: ============================================================================

:: ============================================================================
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
:get_arg_with_value
if "%~3"=="%~2" (
  set "%~1=%~4"
  exit /b 0
)
if not defined %~1 set "%~1=%~5"
exit /b 1

:: ============================================================================
:: Sets a variable from a cli flag argument
:: 1 cli argument name
:: 2 variable name
:: 3 current argument name
:get_arg_flag
if "%~3"=="%~2" (
  set "%~1=1"
  exit /b 0
)
if not defined %~1 set "%~1=0"
exit /b 1

:: ============================================================================
:: https://en.wikipedia.org/wiki/ANSI_escape_code#Select_Graphic_Rendition_parameters
:: https://en.wikipedia.org/wiki/ANSI_escape_code#Colors
:: %ESC%[38;2;<r>;<g>;<b>m TEXT %ESC%[0m
:setESC
for /F "tokens=1,2 delims=#" %%a in ('"prompt #$H#$E# & echo on & for %%b in (1) do rem"') do (
  set ESC=%%b
  exit /B 0
)

:: ============================================================================
:: refs:
:: https://stackoverflow.com/questions/58732724/turn-on-echo-for-a-single-command-in-a-batch-file
:: https://ss64.com/nt/syntax-macros.html
:: https://www.dostips.com/forum/viewtopic.php?f=3&t=2518
:: requires setlocal DisableDelayedExpansion for definition
:define_macro_execute
:: NOTE: Post by Dave Benham (https://www.dostips.com/forum/viewtopic.php?p=58686#p58686)
:: ^^^(<CR>)<LF>(<CR>)<LF>
:: After parsing the definition, the actual stored value is ^<LF>
:: When included at the end of a macro definition line, the ^<LF> is in front of the (<CR>)<LF> at the end of
:: the source line, which is parsed as a single <LF> that gets inserted into the definition, and
:: the next line is still appended to the definition.
(set \n=^^^

)
set ^"execute=for %%# in (1 2) do if %%#==2 ( %\n%
                setlocal EnableDelayedExpansion %\n%
                for /F "tokens=1,*" %%1 in ("!time: =0! !argv!") do ( %\n%
                  @REM echo [%ESC%[32m%%1%ESC%[0m] %%2 %\n%
                  echo [%ESC%[38;2;0;167;204m%%1%ESC%[0m] %%2 ^>^&2%\n%
                  if "%-dry%"=="0" %%2 %\n%
                  endlocal ^& endlocal %\n%
                ) %\n%
              ) else setlocal DisableDelayedExpansion ^& set argv="
exit /b
