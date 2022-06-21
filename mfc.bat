@echo off


if "%1" == ""       goto label_help
if "%1" == "-h"     goto label_help
if "%1" == "--help" goto label_help

goto label_docker


:label_help
echo.
echo mfc.bat: Windows proxy script for mfc.sh using Docker.
echo.
echo   Usage (1): mfc.bat ^[-h^|--help^]      # Print this help message
echo   Usage (2): mfc.bat ^[COMMAND^]        # Equivalent to "./mfc.sh [COMMAND]"
echo   Usage (3): mfc.bat docker           # Obtain an interactive bash session
echo.
echo In all cases, a docker container is used when running a command.
echo.
exit /b 0
:label_help_after


:label_docker
where docker >nul 2>nul
if %errorlevel% neq 0 (
	echo.
	echo ^[mfc.bat^] You must have Docker installed.
	echo ^[mfc.bat^] Please install Docker and try again.
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


if "%1" == "docker" (
	set arg="/bin/bash"
) else (
	set arg="cd MFC && ./mfc.sh %*"
)

echo ^[mfc.bat^] Starting container...
docker run --interactive --tty --rm ^
		   --mount type=bind,source="%cd%",target=/home/me/MFC ^
		   henryleberre/mfc /bin/bash -c %arg%

if %errorlevel% neq 0 (
	echo.
	echo          Docker: Fatal container runtime error.
	exit /b 1
)


exit /b 0
:label_docker_after
