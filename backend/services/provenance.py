"""Build a deployment-neutral provenance and interpretation summary."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


INTERPRETATION_LIMITS = [
    "Scores prioritize hypotheses; they are not experimental validation.",
    "Missing evidence must not be interpreted as evidence of no regulation.",
    "FBA-family predictions depend on model structure, objective, and constraints.",
    "AI-generated recommendations require review against the cited database evidence.",
]


def _table_exists(connection, name: str) -> bool:
    row = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type IN ('table', 'view') AND name = ?",
        (name,),
    ).fetchone()
    return row is not None


def _database_metadata(manager) -> dict[str, Any]:
    connection = manager.get_connection() if manager else None
    if connection is None:
        return {
            "status": "unavailable",
            "schema_version": None,
            "datasets": [],
            "warnings": ["Database metadata is unavailable."],
        }

    schema_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
    datasets = []
    warnings = []
    if _table_exists(connection, "dataset_metadata"):
        rows = connection.execute(
            "SELECT dataset_id, file_path, sha256, byte_size, record_count, "
            "imported_at, source_version, notes FROM dataset_metadata ORDER BY dataset_id"
        ).fetchall()
        datasets = [dict(row) for row in rows]
    else:
        warnings.append(
            "This deployment database omits per-dataset hashes; consult the release manifest."
        )

    migrations = []
    if _table_exists(connection, "schema_migrations"):
        migrations = [
            dict(row)
            for row in connection.execute(
                "SELECT version, description, applied_at FROM schema_migrations ORDER BY version"
            ).fetchall()
        ]

    return {
        "status": "complete" if datasets else "partial",
        "schema_version": schema_version,
        "dataset_count": len(datasets),
        "datasets": datasets,
        "migrations": migrations,
        "warnings": warnings,
    }


def _release_manifest(root: Path) -> dict[str, Any] | None:
    path = root / "data" / "deploy" / "cgl_regulation_public.manifest.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    return {
        "path": path.relative_to(root).as_posix(),
        "bytes": payload.get("bytes"),
        "tables": payload.get("tables"),
        "views": payload.get("views"),
        "row_counts": payload.get("rows", {}),
    }


def _method_versions(root: Path) -> dict[str, str]:
    path = root / "web" / "method_versions.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return {}
    return {
        str(key): str(value)
        for key, value in payload.items()
        if isinstance(key, str) and isinstance(value, (str, int, float))
    }


def build_provenance(*, manager, root: Path, version: str, deployment: str) -> dict[str, Any]:
    """Return auditable metadata without exposing local filesystem paths."""
    database = _database_metadata(manager)
    manifest = _release_manifest(root)
    warnings = list(database.pop("warnings", []))
    if manifest is None:
        warnings.append("The public release manifest is unavailable.")

    return {
        "status": "complete" if database["status"] == "complete" else "partial",
        "app": "cgl-regulation",
        "version": version,
        "deployment": deployment,
        "database": database,
        "release_manifest": manifest,
        "method_versions": _method_versions(root),
        "evidence_policy": {
            "experimental": "Direct or curated experimental evidence.",
            "computed": "Derived by a documented deterministic or statistical workflow.",
            "predicted": "Model-based hypothesis requiring independent validation.",
            "unknown": "Evidence type or provenance is not available.",
        },
        "interpretation_limits": INTERPRETATION_LIMITS,
        "warnings": warnings,
    }
