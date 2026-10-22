@echo off
rem One-time setup for Buddy Builder (Windows). Safe to run again.
rem   setup.bat           install into .venv and check it works
rem   setup.bat --tests   also run the full test suite
cd /d "%~dp0"
echo == Buddy Builder setup ==
where uv >nul 2>nul
if %errorlevel%==0 (
  echo Using uv ^(it downloads Python 3.12 if needed^)
  if not exist .venv\Scripts\python.exe uv venv --python 3.12 .venv || goto :fail
  uv pip install --python .venv\Scripts\python.exe -r requirements.txt || goto :fail
) else (
  if not exist .venv\Scripts\python.exe (
    py -3.12 -m venv .venv 2>nul || python -m venv .venv || goto :nopython
  )
  .venv\Scripts\python.exe -m pip install --quiet --upgrade pip
  .venv\Scripts\python.exe -m pip install -r requirements.txt || goto :fail
)
if not exist .env (
  copy .env.example .env >nul
  echo Created .env: add an AI key there to use an AI model in the chat ^(optional^).
)
echo == Self-check ==
.venv\Scripts\python.exe tools\selfcheck.py || goto :fail
if "%~1"=="--tests" .venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
echo.
echo Ready. Start Buddy Builder with:  start.bat
echo New here? Read GETTING_STARTED.md
exit /b 0
:nopython
echo Python 3.12 was not found. Install it from https://www.python.org/downloads/ (tick "Add python.exe to PATH"),
echo or install uv from https://docs.astral.sh/uv/ , then run setup.bat again.
exit /b 1
:fail
echo Setup failed: see the messages above.
exit /b 1
