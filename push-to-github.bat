@echo off
REM ============================================================
REM  Publish the ICMLSC 2027 (M4338) supplementary materials to
REM  https://github.com/rscheng0707/restart-fair-benchmarking
REM
REM  Put this file in paper2\release\ and double-click it.
REM  Safe to re-run: it skips steps already done.
REM ============================================================
setlocal
cd /d "%~dp0"

echo.
echo === Publishing from: %CD%
echo.

where git >nul 2>&1
if errorlevel 1 (
  echo [X] git is not on PATH. Install Git for Windows, then re-run.
  pause & exit /b 1
)

REM --- sanity: are we in the right folder? ---
if not exist "analysis\reproduce_all.py" (
  echo [X] This does not look like the release folder.
  echo     Expected to find analysis\reproduce_all.py here.
  pause & exit /b 1
)

REM --- git identity, only if not already set ---
for /f "delims=" %%i in ('git config --global user.name 2^>nul') do set GITNAME=%%i
if "%GITNAME%"=="" (
  echo Setting git identity...
  git config --global user.name "Reginald S. Cheng"
  git config --global user.email "rscheng@feutech.edu.ph"
)

REM --- init, only if not already a repo ---
if not exist ".git" (
  echo Initialising repository...
  git init -b main || goto :fail
)

echo Staging 37 files...
git add . || goto :fail

REM --- commit, only if there is something to commit ---
git diff --cached --quiet
if errorlevel 1 (
  echo Committing...
  git commit -F "%~dp0COMMIT_MSG.txt" || goto :fail
) else (
  echo Nothing new to commit.
)

REM --- remote, only if not already set ---
git remote get-url origin >nul 2>&1
if errorlevel 1 (
  echo Adding remote...
  git remote add origin https://github.com/rscheng0707/restart-fair-benchmarking.git || goto :fail
)

echo.
echo Pushing to GitHub. If prompted, sign in with your browser,
echo or paste a personal access token as the password.
echo.
git push -u origin main || goto :fail

echo.
echo ============================================================
echo  DONE. Pushed to:
echo  https://github.com/rscheng0707/restart-fair-benchmarking
echo.
echo  NEXT, in this order -- the order matters:
echo   1. zenodo.org - Account - GitHub - toggle the repo ON
echo   2. THEN on GitHub: Releases - Draft a new release
echo      tag v1.0-icmlsc2027 - Publish release
echo   3. Copy the CONCEPT DOI ("Cite all versions") and send it
echo      to Claude to paste into the paper.
echo ============================================================
echo.
pause & exit /b 0

:fail
echo.
echo [X] A git command failed. The message above says why.
echo     Nothing was lost - you can fix it and re-run this file.
pause & exit /b 1
