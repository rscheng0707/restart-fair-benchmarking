@echo off
REM run_session.bat — runs one main-matrix session with pinning hard-set in the
REM script, not relied on from the calling shell. This exists because session 2
REM ran unpinned: OMP_PROC_BIND/OMP_PLACES were set by hand in the shell and did
REM not carry over. Freeze section 5: pinning is a fixed execution condition,
REM required identical in every session.
REM
REM mainmatrix.exe does not write its own .log file — it prints progress and the
REM pinning confirmation line to stderr only. This script redirects stderr to
REM <out>.log so that record exists to check afterward, same as session1.log /
REM session2.log.
REM
REM Usage:
REM   run_session.bat <session_no> <seeds_from> <seeds_to> <out_csv> [shutdown]
REM
REM   shutdown   optional 5th argument, literal word "shutdown" — powers the
REM              machine off automatically, but ONLY if mainmatrix.exe exits
REM              successfully. A cancelable 60s countdown runs first
REM              (shutdown /a to abort). Omit this argument to just run and
REM              stop, e.g. if you're running sessions back-to-back by hand.
REM
REM Examples:
REM   run_session.bat 2 4 6 session2_pinned.csv
REM   run_session.bat 3 7 9 session3.csv shutdown

setlocal enabledelayedexpansion

if "%~4"=="" (
    echo Usage: run_session.bat ^<session_no^> ^<seeds_from^> ^<seeds_to^> ^<out_csv^> [shutdown]
    exit /b 1
)

set SESSION=%1
set SEEDS_FROM=%2
set SEEDS_TO=%3
set OUT=%4
set DO_SHUTDOWN=%5

set LOGFILE=%OUT:.csv=.log%

REM --- Pinning: fixed condition, not optional, not inherited from the shell ---
set OMP_PROC_BIND=close
set OMP_PLACES=cores

echo ============================================================
echo  session         = %SESSION%
echo  seeds            = %SEEDS_FROM%-%SEEDS_TO%
echo  OMP_PROC_BIND    = %OMP_PROC_BIND%
echo  OMP_PLACES       = %OMP_PLACES%
echo  out              = %OUT%
echo  log              = %LOGFILE%
echo  shutdown after   = %DO_SHUTDOWN%
echo ============================================================

if not "%OMP_PROC_BIND%"=="close" (
    echo FATAL: OMP_PROC_BIND did not get set. Aborting rather than running unpinned.
    exit /b 2
)
if not "%OMP_PLACES%"=="cores" (
    echo FATAL: OMP_PLACES did not get set. Aborting rather than running unpinned.
    exit /b 2
)

mainmatrix.exe --session %SESSION% --seeds-from %SEEDS_FROM% --seeds-to %SEEDS_TO% ^
    --taufile tau.csv --out %OUT% --label main --orderseed 20260916 --warmup 600 ^
    2> "%LOGFILE%"

set EXITCODE=%ERRORLEVEL%

echo.
echo ============================================================
echo  session %SESSION% finished, exit code %EXITCODE%
echo  Log written to %LOGFILE% — confirm it contains:
echo    OMP_PROC_BIND=close OMP_PLACES=cores
echo  before accepting this session as valid.
echo ============================================================

findstr /C:"OMP_PROC_BIND=close OMP_PLACES=cores" "%LOGFILE%" >nul 2>&1
set PIN_OK=%ERRORLEVEL%
if not "%PIN_OK%"=="0" (
    echo WARNING: could not confirm the pinning line in %LOGFILE%.
    echo Do not trust this session, and do NOT shut down automatically.
    endlocal
    exit /b 3
)

if /i "%DO_SHUTDOWN%"=="shutdown" (
    if "%EXITCODE%"=="0" (
        echo.
        echo Session %SESSION% completed successfully and pinning is confirmed in the log.
        echo Shutting down in 60 seconds. Run "shutdown /a" now to cancel.
        shutdown /s /t 60 /c "run_session.bat: session %SESSION% complete, auto-shutdown"
    ) else (
        echo.
        echo Exit code %EXITCODE% is non-zero — NOT shutting down. Check %LOGFILE% first.
    )
)

endlocal
exit /b %EXITCODE%
