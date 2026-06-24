#!/usr/bin/env bash

# Print Header
echo "=========================================================="
echo "   C. glutamicum Regulatory Network Explorer (C.g. RNE)"
echo "                macOS/Linux Startup Launcher"
echo "=========================================================="
echo ""

# 1. Check Python installation
echo "Checking for Python 3..."
if command -v python3 &>/dev/null; then
    PYTHON_CMD="python3"
elif command -v python &>/dev/null; then
    # Fallback to python and check version
    PYTHON_VER=$(python -c 'import sys; print(sys.version_info[0])' 2>/dev/null)
    if [ "$PYTHON_VER" = "3" ]; then
        PYTHON_CMD="python"
    fi
fi

if [ -z "$PYTHON_CMD" ]; then
    echo "[ERROR] Python 3 was not found on your system!"
    echo ""
    echo "Please install Python 3 using your package manager (brew, apt, pacman) or download it from:"
    echo "https://www.python.org/downloads/"
    echo ""
    read -p "Press Enter to exit..."
    exit 1
fi

PY_VERSION_STR=$($PYTHON_CMD --version)
echo "Found $PY_VERSION_STR"

# 2. Install dependencies
echo ""
echo "Installing/upgrading required Python libraries..."
$PYTHON_CMD -m pip install --upgrade pip
$PYTHON_CMD -m pip install -r requirements.txt

if [ $? -ne 0 ]; then
    echo ""
    echo "[WARNING] Failed to install dependencies automatically."
    echo "If this is your first time running, some features might not work properly."
    echo "You may need to run: pip3 install pandas networkx matplotlib pyvis"
    echo ""
    read -p "Do you want to try starting the server anyway? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# 3. Start local web server
echo ""
echo "Starting local web server on port 8000..."
echo "The network explorer will open in your default browser shortly."
echo "To stop the server, press Ctrl+C in this terminal."
echo "----------------------------------------------------------"
$PYTHON_CMD run_server.py

if [ $? -ne 0 ]; then
    echo ""
    echo "[ERROR] Local server crashed or failed to start."
    read -p "Press Enter to exit..."
fi
