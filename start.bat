@echo off
rem Start Buddy Builder as an app in the browser (Windows). Options: see start.py.
cd /d "%~dp0"
if exist .venv\Scripts\python.exe (
  .venv\Scripts\python.exe start.py %*
) else (
  python start.py %*
)
