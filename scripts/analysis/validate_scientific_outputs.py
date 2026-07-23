#!/usr/bin/env python3
"""Generate a machine-readable audit of scientific output integrity.

This checks traceability and numerical invariants. It deliberately does not
claim biological validation, which requires an independent held-out dataset.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from pathlib import Path


SCORE_COLUMNS = {
    "regulations": ["evidence_score"],
    "network_edges_extended": ["confidence_score"],
    "condition_edge_scores": [
        "base_confidence", "dynamic_support", "evidence_completeness", "condition_score"
    ],
    "intervention_target_scores": [
        "evidence_score", "systems_impact_score", "engineering_tractability_score",
        "risk_score", "priority_score",
    ],
    "tf_gene_rf_scores": ["predicted_confidence", "evidence_score"],
}


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average = ((start + 1) + end) / 2
        for position in range(start, end):
            ranks[order[position]] = average
        start = end
    return ranks


def evaluate_gold_standard(
    path: Path, *, score_column: str = "score", label_column: str = "label", bins: int = 10
) -> dict:
    """Evaluate binary predictions without using the training-time ML stack."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("Gold-standard CSV contains no rows")
    required = {score_column, label_column}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Gold-standard CSV is missing columns: {', '.join(sorted(missing))}")

    scores = []
    labels = []
    for line_number, row in enumerate(rows, start=2):
        try:
            score = float(row[score_column])
            label = int(row[label_column])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid score or label at CSV line {line_number}") from exc
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError(f"Score outside [0, 1] at CSV line {line_number}")
        if label not in {0, 1}:
            raise ValueError(f"Label must be 0 or 1 at CSV line {line_number}")
        scores.append(score)
        labels.append(label)

    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise ValueError("Gold-standard CSV must contain both positive and negative labels")

    ranks = _average_ranks(scores)
    positive_rank_sum = sum(rank for rank, label in zip(ranks, labels) if label == 1)
    auroc = (positive_rank_sum - positives * (positives + 1) / 2) / (positives * negatives)

    ordered = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    true_positives = 0
    precision_sum = 0.0
    for rank, (_score, label) in enumerate(ordered, start=1):
        if label:
            true_positives += 1
            precision_sum += true_positives / rank
    average_precision = precision_sum / positives
    brier_score = sum((score - label) ** 2 for score, label in zip(scores, labels)) / len(labels)

    calibration = []
    expected_calibration_error = 0.0
    for index in range(bins):
        lower = index / bins
        upper = (index + 1) / bins
        members = []
        for score, label in zip(scores, labels):
            in_bin = lower <= score <= upper if index == bins - 1 else lower <= score < upper
            if in_bin:
                members.append((score, label))
        if not members:
            continue
        mean_score = sum(item[0] for item in members) / len(members)
        observed_rate = sum(item[1] for item in members) / len(members)
        expected_calibration_error += len(members) / len(labels) * abs(mean_score - observed_rate)
        calibration.append({
            "lower": lower, "upper": upper, "count": len(members),
            "mean_score": mean_score, "observed_rate": observed_rate,
        })

    return {
        "status": "evaluated",
        "file": path.name,
        "score_column": score_column,
        "label_column": label_column,
        "samples": len(labels),
        "positives": positives,
        "negatives": negatives,
        "metrics": {
            "auroc": round(auroc, 6),
            "average_precision": round(average_precision, 6),
            "brier_score": round(brier_score, 6),
            "expected_calibration_error": round(expected_calibration_error, 6),
        },
        "calibration_bins": calibration,
    }


def _table_exists(connection: sqlite3.Connection, name: str) -> bool:
    return connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (name,)
    ).fetchone() is not None


def audit(
    database: Path, *, gold_standard: Path | None = None,
    score_column: str = "score", label_column: str = "label",
) -> dict:
    connection = sqlite3.connect(f"file:{database.resolve().as_posix()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    checks = []
    try:
        schema_version = int(connection.execute("PRAGMA user_version").fetchone()[0])
        checks.append({"name": "schema_version", "passed": schema_version > 0, "value": schema_version})

        for table, columns in SCORE_COLUMNS.items():
            if not _table_exists(connection, table):
                checks.append({"name": f"{table}.present", "passed": False})
                continue
            count = int(connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
            checks.append({"name": f"{table}.nonempty", "passed": count > 0, "value": count})
            available = {row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')}
            for column in columns:
                if column not in available:
                    checks.append({"name": f"{table}.{column}.present", "passed": False})
                    continue
                row = connection.execute(
                    f'SELECT MIN("{column}"), MAX("{column}"), '
                    f'SUM(CASE WHEN "{column}" < 0 OR "{column}" > 1 THEN 1 ELSE 0 END) '
                    f'FROM "{table}" WHERE "{column}" IS NOT NULL'
                ).fetchone()
                checks.append({
                    "name": f"{table}.{column}.unit_interval",
                    "passed": int(row[2] or 0) == 0,
                    "min": row[0], "max": row[1], "out_of_range": int(row[2] or 0),
                })

        if _table_exists(connection, "regulations"):
            total, with_pmid = connection.execute(
                "SELECT COUNT(*), SUM(CASE WHEN PMID IS NOT NULL AND TRIM(PMID) != '' THEN 1 ELSE 0 END) "
                "FROM regulations"
            ).fetchone()
            checks.append({
                "name": "regulations.pmid_coverage",
                "passed": total > 0,
                "value": round((with_pmid or 0) / total, 4) if total else 0,
                "informational": True,
            })

        dataset_count = 0
        if _table_exists(connection, "dataset_metadata"):
            dataset_count = int(connection.execute("SELECT COUNT(*) FROM dataset_metadata").fetchone()[0])
            invalid_hashes = int(connection.execute(
                "SELECT COUNT(*) FROM dataset_metadata WHERE LENGTH(sha256) != 64"
            ).fetchone()[0])
            checks.append({
                "name": "dataset_metadata.sha256_coverage",
                "passed": dataset_count > 0 and invalid_hashes == 0,
                "datasets": dataset_count, "invalid_hashes": invalid_hashes,
            })
    finally:
        connection.close()

    failed = [item["name"] for item in checks if not item["passed"] and not item.get("informational")]
    report = {
        "status": "pass" if not failed else "fail",
        "database": database.name,
        "checks": checks,
        "failed_checks": failed,
        "scope": "structural integrity, traceability, and score-domain validation",
        "not_validated": [
            "biological accuracy against an independent gold-standard dataset",
            "causal effects of proposed interventions",
            "wet-lab efficacy of AI-generated recommendations",
        ],
    }
    if gold_standard:
        report["external_validation"] = evaluate_gold_standard(
            gold_standard, score_column=score_column, label_column=label_column
        )
        report["external_validation"]["interpretation"] = (
            "Metrics are valid only if the supplied rows are independent of training and feature construction."
        )
    else:
        report["external_validation"] = {
            "status": "not_provided",
            "required_columns": [label_column, score_column],
        }
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, default=Path("data/reference/cgl_regulation.db"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--gold-standard", type=Path)
    parser.add_argument("--score-column", default="score")
    parser.add_argument("--label-column", default="label")
    args = parser.parse_args()
    try:
        report = audit(
            args.database,
            gold_standard=args.gold_standard,
            score_column=args.score_column,
            label_column=args.label_column,
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    rendered = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    raise SystemExit(0 if report["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
