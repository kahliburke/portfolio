@echo off
rem Double-click launcher for run.ps1. PowerShell holds the launch logic; this only
rem exists so a double-click runs it (a .ps1 double-click opens an editor, not runs).
set SLATE_LAUNCHED_BY_BAT=1
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run.ps1"
rem Backstop: keep the window open even if PowerShell itself was blocked from running the
rem script (e.g. locked-down machine), so the error is readable rather than flashing shut.
pause
