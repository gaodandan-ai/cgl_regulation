@echo off
title C. glutamicum Regulatory Network Explorer
echo ==========================================================
echo    C. glutamicum Regulatory Network Explorer (C.g. RNE)
echo                  Windows Startup Launcher
echo ==========================================================
echo.

:: 1. Check Python installation and find the correct command
echo Checking for Python 3...

:: Test 'py -3' (Python Launcher)
py -3 -c "import sys" >nul 2>nul
if %errorlevel% equ 0 (
    set PY_CMD=py -3
    goto :found_python
)

:: Test standard 'python' command (checking if it runs properly, not just 'where')
python -c "import sys; sys.exit(0 if sys.version_info[0] == 3 else 1)" >nul 2>nul
if %errorlevel% equ 0 (
    set PY_CMD=python
    goto :found_python
)

:: If we reach here, Python was not found or is misconfigured
echo [ERROR] Python 3 was not found or is not properly configured on your system.
echo.
echo Please do one of the following:
echo 1. Download and install Python 3.x from:
echo    https://www.python.org/downloads/
echo    (Make sure to check "Add Python to PATH" during installation)
echo.
echo 2. If already installed, disable the Microsoft Store Python App Alias:
echo    Go to: Settings > Apps > Advanced app settings > App execution aliases
echo    And turn OFF the switch for "Python" and "Python3".
echo.
pause
exit /b 1

:found_python
:: Show Python version
for /f "tokens=*" %%i in ('%PY_CMD% --version') do set py_ver=%%i
echo Found %py_ver% (Using command: %PY_CMD%)

:: 2. Install dependencies
echo.
echo Installing/upgrading required Python libraries...
%PY_CMD% -m pip install --upgrade pip
%PY_CMD% -m pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo.
    echo [WARNING] Failed to install some dependencies automatically.
    echo If this is your first time running, some features might not work properly.
    echo You may need to run: pip install pandas networkx matplotlib pyvis
    echo.
    choice /c yn /m "Do you want to try starting the server anyway?"
    if errorlevel 2 (
        exit /b 1
    )
)

:: 3. Start local web server
echo.
echo Starting local web server on port 8000...
echo The network explorer will open in your default browser shortly.
echo To stop the server, close this window or press Ctrl+C in this terminal.
echo ----------------------------------------------------------
%PY_CMD% run_server.py

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] Local server crashed or failed to start.
    pause
)
