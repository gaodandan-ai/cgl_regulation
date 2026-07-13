"""
scripts/build_app.py
====================
Packages the Cgl Regulation explorer platform into a single standalone
desktop executable (dist/cgl_regulation.exe) using PyInstaller.

Includes all static files, models, and hidden imports for uvicorn & scikit-learn.
"""

import os
import subprocess
import sys

def main():
    # Go to repository root
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(ROOT)

    print("=== Preparing PyInstaller build ===")

    # Setup PyInstaller arguments
    cmd = [
        "pyinstaller",
        "--onefile",
        "--name", "cgl_regulation",
        "--icon", "icon.ico",
        # Add backend/ to paths during analysis
        "--paths", "backend",
        # Add static assets and data folders
        "--add-data", "web;web",
        "--add-data", "data/reference;data/reference",
        "--add-data", "backend/models;backend/models",
        # Uvicorn hidden imports
        "--hidden-import", "uvicorn.protocols.http.h11_impl",
        "--hidden-import", "uvicorn.loop.asyncio",
        "--hidden-import", "uvicorn.lifespan.on",
        "--hidden-import", "uvicorn.lifespan.off",
        "--hidden-import", "uvicorn.loop.auto",
        # Scikit-learn forest models & dependencies
        "--hidden-import", "sklearn.ensemble._forest",
        "--hidden-import", "sklearn.utils._typedefs",
        "--hidden-import", "sklearn.neighbors._typedefs",
        # FastAPI / Cobra / Depinfo
        "--hidden-import", "fastapi",
        "--hidden-import", "pydantic",
        "--hidden-import", "cobra",
        "--hidden-import", "depinfo",
        # RAG service
        "--hidden-import", "rag_service",
        # Backend modules (both absolute and relative names due to sys.path manipulation)
        "--hidden-import", "backend.app",
        "--hidden-import", "backend.gene_utils",
        "--hidden-import", "backend.kegg_client",
        "--hidden-import", "backend.metabolic_mapper",
        "--hidden-import", "backend.bio_handlers",
        "--hidden-import", "backend.sequence_tools",
        "--hidden-import", "backend.model_loader",
        "--hidden-import", "backend.thermo_pruner",
        "--hidden-import", "backend.simulation",
        "--hidden-import", "backend.schemas",
        "--hidden-import", "backend.objectives",
        "--hidden-import", "backend.thermodynamics",
        "--hidden-import", "backend.enzyme_thermal_params",
        "--hidden-import", "app",
        "--hidden-import", "gene_utils",
        "--hidden-import", "kegg_client",
        "--hidden-import", "metabolic_mapper",
        "--hidden-import", "bio_handlers",
        "--hidden-import", "sequence_tools",
        "--hidden-import", "model_loader",
        "--hidden-import", "thermo_pruner",
        "--hidden-import", "simulation",
        "--hidden-import", "schemas",
        "--hidden-import", "objectives",
        "--hidden-import", "thermodynamics",
        "--hidden-import", "enzyme_thermal_params",
        # Main entrypoint script
        "run_server.py"
    ]

    print(f"Running command: {' '.join(cmd)}\n")
    try:
        subprocess.run(cmd, check=True)
        print("\n=== PyInstaller build completed successfully ===")
        print("Executable generated at: dist/cgl_regulation.exe")
    except subprocess.CalledProcessError as e:
        print(f"\nError: PyInstaller build failed with exit code {e.returncode}")
        sys.exit(e.returncode)

if __name__ == '__main__':
    main()
