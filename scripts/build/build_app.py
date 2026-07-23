#!/usr/bin/env python3
"""
scripts/build_app.py
====================
CI build script for GitHub Actions.
- Reads the current git tag as the version string
- Updates web/version.json with the new version
- Runs PyInstaller to produce dist/cgl_regulation.exe
"""

import json
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def get_version():
    """Return version from RELEASE_VERSION env var (set by GitHub Actions), git tag, or fallback."""
    # GitHub Actions sets this from the tag
    env_ver = os.environ.get("RELEASE_VERSION", "").strip().lstrip("v")
    if env_ver:
        return env_ver

    # Try git describe
    try:
        tag = subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0"],
            cwd=ROOT, stderr=subprocess.DEVNULL
        ).decode().strip().lstrip("v")
        if tag:
            return tag
    except Exception:
        pass

    return "0.0.0"


def update_version_json(version: str):
    """Write version info to web/version.json."""
    import datetime
    path = os.path.join(ROOT, "web", "version.json")
    data = {
        "version": version,
        "release_date": datetime.date.today().isoformat(),
        "download_url": "https://github.com/gaodandan-ai/cgl_regulation/releases/latest/download/cgl_regulation.exe",
        "changelog": os.environ.get("RELEASE_NOTES", "See GitHub release notes for details.")
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"[build] Updated web/version.json -> v{version}")


def write_local_version(version: str):
    """Write a bundled version file that the frozen exe can read at runtime."""
    path = os.path.join(ROOT, "web", "version_local.json")
    data = {"version": version}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
    print(f"[build] Wrote web/version_local.json -> v{version}")


def run_pyinstaller():
    spec = os.path.join(ROOT, "cgl_regulation.spec")
    result = subprocess.run(
        [sys.executable, "-m", "PyInstaller", spec, "--noconfirm"],
        cwd=ROOT
    )
    if result.returncode != 0:
        print("[build] PyInstaller failed!", file=sys.stderr)
        sys.exit(1)
    print("[build] PyInstaller build succeeded.")


if __name__ == "__main__":
    version = get_version()
    print(f"[build] Building version v{version}")
    update_version_json(version)
    write_local_version(version)
    run_pyinstaller()
    print(f"[build] Done. Output: dist/cgl_regulation.exe")
