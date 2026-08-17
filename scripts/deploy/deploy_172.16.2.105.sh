#!/usr/bin/env bash
# ==============================================================================
# deploy_172.16.2.105.sh
# ==============================================================================
# One-click Intranet Server Deployment & Update Script for 172.16.2.105
# C. glutamicum Regulatory Network Explorer (Cgl Regulation Explorer)
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SERVER_IP="172.16.2.105"
PORT=8010
COMPOSE_FILE="${SCRIPT_DIR}/docker-compose.intranet.yml"

echo "======================================================================"
echo "   C. glutamicum Regulatory Network Explorer (Cgl Regulation Explorer)"
echo "            Intranet Server One-Click Deployment (${SERVER_IP})"
echo "======================================================================"
echo ""

# Function: check docker availability
check_docker() {
    if command -v docker &>/dev/null && docker info &>/dev/null; then
        return 0
    else
        return 1
    fi
}

# 1. Option parsing
ACTION="${1:-deploy}"

case "$ACTION" in
    status)
        echo "Checking service status on ${SERVER_IP}:${PORT}..."
        if check_docker; then
            docker compose -f "$COMPOSE_FILE" ps
        fi
        curl -s "http://127.0.0.1:${PORT}/api/health" || echo "[OFFLINE] Service is not responding on port ${PORT}"
        exit 0
        ;;
    stop)
        echo "Stopping intranet service..."
        if check_docker; then
            docker compose -f "$COMPOSE_FILE" down
        else
            pkill -f "run_server.py" || true
        fi
        echo "Service stopped."
        exit 0
        ;;
    logs)
        echo "Viewing server logs..."
        if check_docker; then
            docker compose -f "$COMPOSE_FILE" logs -f --tail=100
        else
            tail -n 100 ~/cgl_server.log || true
        fi
        exit 0
        ;;
esac

# 2. Main deployment logic
echo "[1/4] Preparing data and environment..."
python3 data_pipeline/scripts/import_lab_chip_seq_edges.py || true
python3 data_pipeline/scripts/import_lab_expression_compendium.py || true
python3 data_pipeline/scripts/import_lab_chip_peaks.py || true
python3 data_pipeline/scripts/build_sqlite_db.py || true

if check_docker; then
    echo "[2/4] Deploying with Docker Compose (${COMPOSE_FILE})..."
    docker compose -f "$COMPOSE_FILE" up -d --build
    echo "[3/4] Waiting for container health check..."
    sleep 5
    docker compose -f "$COMPOSE_FILE" ps
else
    echo "[WARNING] Docker not found or permission denied. Falling back to native Python..."
    echo "[2/4] Installing dependencies..."
    python3 -m pip install -r requirements-core.txt -r requirements.txt || true
    echo "[3/4] Starting server in background..."
    export PORT=8010
    export CGL_HOST=0.0.0.0
    export HEADLESS=true
    export CGL_PUBLIC_DEPLOYMENT=false
    nohup python3 run_server.py > ~/cgl_server.log 2>&1 &
fi

echo "[4/4] Verifying deployment on http://${SERVER_IP}:${PORT}..."
sleep 3
if curl -s "http://127.0.0.1:${PORT}/api/health" | grep -q "cgl-regulation"; then
    echo ""
    echo "======================================================================"
    echo " SUCCESS! Intranet Server is successfully running at:"
    echo " 👉 http://${SERVER_IP}:${PORT}"
    echo " 👉 http://${SERVER_IP}:${PORT}/index.html"
    echo "======================================================================"
else
    echo ""
    echo "[WARNING] Server started, but health check endpoint returned non-200."
    echo "Check logs using: bash deploy_172.16.2.105.sh logs"
fi
