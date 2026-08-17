@echo off
REM ==============================================================================
REM deploy_172.16.2.105.bat
REM ==============================================================================
REM One-click Intranet Windows Server Deployment & Startup Script for 172.16.2.105
REM C. glutamicum Regulatory Network Explorer (Cgl Regulation Explorer)
REM ==============================================================================

echo ======================================================================
echo    C. glutamicum Regulatory Network Explorer (Cgl Regulation Explorer)
echo             Intranet Windows Server Launcher (172.16.2.105)
echo ======================================================================
echo.

cd /d "%~dp0..\.."

set PORT=8010
set CGL_HOST=0.0.0.0
set HEADLESS=true
set CGL_PUBLIC_DEPLOYMENT=false
set CGL_INTRANET_SERVER=172.16.2.105

echo [1/3] Consolidating internal omics datasets...
python data_pipeline/scripts/import_lab_chip_seq_edges.py
python data_pipeline/scripts/import_lab_expression_compendium.py
python data_pipeline/scripts/import_lab_chip_peaks.py
python data_pipeline/scripts/build_sqlite_db.py

echo.
echo [2/3] Checking Docker availability...
docker info >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo Docker detected! Starting intranet service container...
    docker compose -f scripts/deploy/docker-compose.intranet.yml up -d --build
    echo Service successfully launched in Docker!
    echo Access at: http://172.16.2.105:8010
    pause
    exit /b 0
)

echo [3/3] Docker not running. Starting native Python Uvicorn server on 0.0.0.0:8010...
echo.
echo ======================================================================
echo  Server running at: http://172.16.2.105:8010/index.html
echo  To stop the server, close this window.
echo ======================================================================
echo.

python run_server.py --port 8010 --no-browser
pause
