@echo off
setlocal EnableDelayedExpansion

call "%~dp0get_conda_exec.bat"
if !errorlevel! neq 0 (
  exit /B !errorlevel!
)

set MY_CONDA=!MY_CONDA_EXE:"=!
cd %~dp0
set PYTHONUTF8=1
call "!MY_CONDA!" remove --name gams --all --yes
call "!MY_CONDA!" env update -f environment-3.9.yml -n gams
call "!MY_CONDA!" activate gams
call python -m pip install -e . --no-deps
pause
cmd /k
