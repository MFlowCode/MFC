@echo off


if "%1" == ""       goto label_help
if "%1" == "-h"     goto label_help
if "%1" == "--help" goto label_help
if "%1" == "docker" goto label_docker

goto label_windows


:label_help
echo.
echo mfc.bat: Windows proxy script for mfc.sh using Docker.
echo.
echo   Usage (1): mfc.bat ^[COMMAND^]        # Equivalent to "./mfc.sh [COMMAND]"
echo   Usage (2): mfc.bat docker           # Obtain an interactive bash session
echo.
exit /b 0
:label_help_after


:label_windows

if not exist "%cd%\toolchain\main.py" (
  echo.
  echo ^[mfc.bat^] You must call this script from within MFC's root folder
  echo.
  exit /b %errorlevel%
)

mkdir "%cd%\build" 2> NUL

if not exist "%cd%\build\venv" (
	python3 -m venv "%cd%\build\venv"
	if %errorlevel% neq 0 (
		echo.
		echo ^[mfc.bat^] Failed to create the Python virtual environment. Delete the build/venv folder and try again.
		echo.
		exit /b %errorlevel%
	)

	pip3 install pyyaml fypp rich argparse dataclasses 
	if %errorlevel% neq 0 (
		echo.
		echo ^[mfc.bat^] Failed to install Python dependencies via Pip. Delete the build/venv folder and try again.
		echo.
		exit /b %errorlevel%
	)
)

call "%cd%\build\venv\Scripts\activate.bat"
if %errorlevel% neq 0 (
	echo.
	echo ^[mfc.bat^] Failed to activate the Python virtual environment.
	echo.
	exit /b %errorlevel%
)

python3 "%cd%\toolchain\main.py" %*
set main_py_err=%errorlevel%

call "%cd%\build\venv\Scripts\deactivate.bat"

if %main_py_err% neq 0 (
	echo.
	echo ^[mfc.bat^] Failed to run MFC.
	echo.
)

exit /b %main_py_err%

:label_windows_after


:label_docker
where docker >nul 2>nul
if %errorlevel% neq 0 (
	echo.
	echo ^[mfc.bat^] You must have Docker installed.
	echo           Please install Docker and try again.
	exit /b %errorlevel%
)


echo ^[mfc.bat^] Fetching image...
docker pull henryleberre/mfc
if %errorlevel% neq 0 (
	echo.
	echo ^[mfc.bat^] Docker: Failed to fetch image.
	echo           Pleasure ensure docker is running.
	exit /b %errorlevel%
)

echo ^[mfc.bat^] Starting container...
docker run --interactive --tty --rm ^
		   --mount type=bind,source="%cd%",target=/home/me/MFC ^
		   henryleberre/mfc

if %errorlevel% neq 0 (
	echo.
	echo          Docker: Fatal container runtime error.
	exit /b %errorlevel%
)


exit /b 0
:label_docker_after
