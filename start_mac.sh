#!/bin/bash
# ====================================================================
# Cgl Regulation Explorer — macOS Launch Script (v1.3.0)
# ====================================================================

echo "============================================================"
echo "  Cgl Regulation Explorer — Starting macOS Desktop Client"
echo "============================================================"

# Navigate to project directory
CDPATH="" cd -- "$(dirname -- "$0")"

# Check Python 3
if command -v python3 &>/dev/null; then
    PYTHON_CMD="python3"
elif command -v python &>/dev/null; then
    PYTHON_CMD="python"
else
    echo "[Error] Python 3 is not installed or not in PATH."
    echo "Please install Python 3.10+ via Homebrew or from https://python.org"
    exit 1
fi

echo "[Info] Using Python interpreter: $($PYTHON_CMD --version)"

# Launch Python Desktop Launcher
$PYTHON_CMD launcher.pyw
