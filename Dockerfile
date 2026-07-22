# ─── Stage 1: Dependencies ────────────────────────────────────────────────────
FROM python:3.11-slim AS deps

WORKDIR /app
COPY requirements.txt .

# Install runtime dependencies (exclude dev/build-only packages)
# equilibrator-api is only needed to rebuild thermo_dgr_data.json offline;
# the pre-built JSON is already in the repo, so we skip it for the runtime image.
RUN pip install --no-cache-dir \
    pandas openpyxl networkx matplotlib pyvis \
    scikit-learn joblib \
    "cobra>=0.22.0" "fastapi>=0.95.0" "uvicorn>=0.20.0" "pydantic>=2.0"

# ─── Stage 2: Runtime image ───────────────────────────────────────────────────
FROM python:3.11-slim

WORKDIR /app

# Copy installed packages from deps stage
COPY --from=deps /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=deps /usr/local/bin /usr/local/bin

# Copy application code and data
COPY . .

# Expose port 8000
EXPOSE 8000

# Environment variables
# PORT: allow hosting platforms (Render, Heroku) to override
# HEADLESS: disable browser auto-open in headless cloud deployments
ENV PORT=8000
ENV HEADLESS=true
ENV CGL_HOST=0.0.0.0
ENV CGL_PUBLIC_DEPLOYMENT=true

# Health check (optional but good practice for orchestration)
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:${PORT}/api/model/status')" || exit 1

# Command to run the application
CMD ["python", "run_server.py"]
