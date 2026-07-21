@echo off
REM 启动客户端 — 无控制台窗口，调用 launcher.pyw
REM 使用 pythonw.exe (Windows 无窗口 Python) 启动，不显示任何黑色窗口
start "" /B pythonw "%~dp0launcher.pyw"
