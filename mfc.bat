@echo off

if "%1" == "docker" goto label_docker

goto label_windows

:label_windows

if not exist "%cd%\toolchain\mfc.py" (
  echo.
  echo ^[mfc.bat^] You must call this script from within MFC's root folder
  echo.
  exit /b 1
)

mkdir "%cd%\build" 2> NUL

if not exist "%cd%\build\venv" (
	python3 -m venv "%cd%\build\venv"
	if %errorlevel% neq 0 (
		echo.
		echo ^[mfc.bat^] Failed to create the Python virtual environment. Delete the build/venv folder and try again.
		echo.
		exit /b 1
	)
)

call "%cd%\build\venv\Scripts\activate.bat"
if %errorlevel% neq 0 (
	echo.
	echo ^[mfc.bat^] Failed to activate the Python virtual environment.
	echo.
	exit /b 1
)

fc /b "%cd%\build\requirements.txt" "%cd%\toolchain\requirements.txt" 2> NUL
if %errorlevel% neq 0 (
    pip3 install -r toolchain/requirements.txt

    copy "%cd%\toolchain\requirements.txt" "%cd%\build" 2> NUL
)

python3 "%cd%\toolchain\mfc.py" %*
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
	exit /b 1
)


echo ^[mfc.bat^] Fetching image...
docker pull henryleberre/mfc
if %errorlevel% neq 0 (
	echo.
	echo ^[mfc.bat^] Docker: Failed to fetch image.
	echo           Pleasure ensure docker is running.
	exit /b 1
)

echo ^[mfc.bat^] Starting container...
docker run --interactive --tty --rm ^
		   --mount type=bind,source="%cd%",target=/home/me/MFC ^
		   henryleberre/mfc

if %errorlevel% neq 0 (
	echo.
	echo          Docker: Fatal container runtime error.
	exit /b 1
)


exit /b 0
:label_docker_after
