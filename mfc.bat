@echo off

goto label_windows

:label_windows

if not exist "%cd%\toolchain\main.py" (
  echo.
  echo ^[mfc.bat^] You must call this script from within MFC's root folder
  echo.
  exit /b 1
)

mkdir "%cd%\build" 2> NUL

if not exist "%cd%\build\venv" (
	python -m venv "%cd%\build\venv"
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

fc /b "%cd%\build\pyproject.toml" "%cd%\toolchain\pyproject.toml" 2> NUL
if %errorlevel% neq 0 (
    pip3 install -e toolchain

    copy "%cd%\toolchain\pyproject.toml" "%cd%\build" 2> NUL
)

python "%cd%\toolchain\main.py" %*
set main_py_err=%errorlevel%

call "%cd%\build\venv\Scripts\deactivate.bat"

if %main_py_err% neq 0 (
	echo.
	echo ^[mfc.bat^] Failed to run MFC.
	echo.
)

exit /b %main_py_err%

:label_windows_after
