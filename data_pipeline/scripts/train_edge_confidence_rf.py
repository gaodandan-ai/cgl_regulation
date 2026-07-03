#!/usr/bin/env python3
"""Train a random-forest TF-target edge confidence model.

Input:
    data/edge_confidence/tf_gene_edge_features.csv

Outputs:
    data/edge_confidence/tf_gene_edge_scores.csv
    data/edge_confidence/edge_confidence_metrics.json
    data/edge_confidence/edge_confidence_feature_importance.csv
    data/edge_confidence/edge_confidence_rf.joblib

By default the script excludes direct curated-evidence fields such as Evidence,
PMID, Source, and regulation role to reduce label leakage. Use
--include-direct-evidence only for exploratory evidence-completeness scoring.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "edge_confidence" / "tf_gene_edge_features.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data" / "edge_confidence"

ID_COLUMNS = {
    "tf_locus",
    "tf_name",
    "target_locus",
    "target_name",
    "label",
    "sample_type",
}

TEXT_COLUMNS = {
    "evidence_text",
    "source",
    "target_operon",
}

DEFAULT_LEAKAGE_COLUMNS = {
    "regulation_role",
    "is_activation",
    "is_repression",
    "is_sigma_factor",
    "evidence_score",
    "has_experimental_evidence",
    "has_pmid",
    "pmid_count",
    "has_database_source",
    "source_coryneregnet",
}


def require_ml_dependencies():
    try:
        import joblib  # noqa: F401
        import sklearn  # noqa: F401
    except Exception as exc:
        raise SystemExit(
            "Missing machine-learning dependencies. Install them with:\n"
            "  python -m pip install scikit-learn joblib\n"
            "or run:\n"
            "  python -m pip install -r requirements.txt\n"
            f"\nOriginal import error: {exc}"
        )


def safe_float(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except Exception:
        return None


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def choose_feature_columns(rows: Sequence[Dict[str, str]], include_direct_evidence: bool) -> List[str]:
    if not rows:
        return []
    excluded = set(ID_COLUMNS) | set(TEXT_COLUMNS)
    if not include_direct_evidence:
        excluded |= DEFAULT_LEAKAGE_COLUMNS

    columns = []
    for column in rows[0].keys():
        if column in excluded:
            continue
        values = [safe_float(row.get(column)) for row in rows]
        numeric_values = [value for value in values if value is not None]
        if numeric_values:
            columns.append(column)
    return columns


def matrix_from_rows(rows: Sequence[Dict[str, str]], feature_columns: Sequence[str]) -> List[List[float]]:
    matrix: List[List[float]] = []
    for row in rows:
        values = []
        for column in feature_columns:
            value = safe_float(row.get(column))
            values.append(value if value is not None else float("nan"))
        matrix.append(values)
    return matrix


def labels_from_rows(rows: Sequence[Dict[str, str]]) -> List[int]:
    labels = []
    for row in rows:
        parsed = safe_float(row.get("label"))
        labels.append(1 if parsed and parsed > 0 else 0)
    return labels


def split_indices(rows: Sequence[Dict[str, str]], labels: Sequence[int], test_size: float, seed: int) -> Tuple[List[int], List[int], str]:
    from sklearn.model_selection import GroupShuffleSplit, train_test_split

    indices = list(range(len(rows)))
    groups = [row.get("tf_locus", "") or f"row_{i}" for i, row in enumerate(rows)]
    unique_groups = sorted(set(groups))
    if len(unique_groups) >= 4:
        splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
        train_idx, test_idx = next(splitter.split(indices, labels, groups))
        return list(train_idx), list(test_idx), "GroupShuffleSplit by tf_locus"

    train_idx, test_idx = train_test_split(
        indices,
        test_size=test_size,
        random_state=seed,
        stratify=labels if len(set(labels)) == 2 else None,
    )
    return list(train_idx), list(test_idx), "train_test_split fallback"


def metric_or_none(fn, *args, **kwargs):
    try:
        result = fn(*args, **kwargs)
        return float(result) if result is not None else None
    except Exception:
        return None


def train_model(args: argparse.Namespace) -> Dict[str, Any]:
    require_ml_dependencies()

    import joblib
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.metrics import (
        accuracy_score,
        average_precision_score,
        brier_score_loss,
        f1_score,
        precision_score,
        recall_score,
        roc_auc_score,
    )
    from sklearn.pipeline import Pipeline

    rows = read_rows(args.input)
    if not rows:
        raise SystemExit(f"No rows found in {args.input}")

    feature_columns = choose_feature_columns(rows, args.include_direct_evidence)
    if not feature_columns:
        raise SystemExit("No numeric feature columns found.")

    x = matrix_from_rows(rows, feature_columns)
    y = labels_from_rows(rows)
    train_idx, test_idx, split_strategy = split_indices(rows, y, args.test_size, args.seed)

    x_train = [x[i] for i in train_idx]
    y_train = [y[i] for i in train_idx]
    x_test = [x[i] for i in test_idx]
    y_test = [y[i] for i in test_idx]

    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("rf", RandomForestClassifier(
            n_estimators=args.n_estimators,
            min_samples_leaf=args.min_samples_leaf,
            class_weight="balanced_subsample",
            random_state=args.seed,
            n_jobs=args.n_jobs,
        )),
    ])
    pipeline.fit(x_train, y_train)

    test_scores = pipeline.predict_proba(x_test)[:, 1].tolist()
    test_pred = [1 if score >= args.threshold else 0 for score in test_scores]
    all_scores = pipeline.predict_proba(x)[:, 1].tolist()

    metrics = {
        "input": str(args.input),
        "row_count": len(rows),
        "positive_count": int(sum(y)),
        "negative_count": int(len(y) - sum(y)),
        "feature_count": len(feature_columns),
        "feature_columns": feature_columns,
        "excluded_columns": sorted((set(rows[0].keys()) - set(feature_columns))),
        "include_direct_evidence": bool(args.include_direct_evidence),
        "split_strategy": split_strategy,
        "train_rows": len(train_idx),
        "test_rows": len(test_idx),
        "threshold": args.threshold,
        "roc_auc": metric_or_none(roc_auc_score, y_test, test_scores),
        "average_precision": metric_or_none(average_precision_score, y_test, test_scores),
        "brier_score": metric_or_none(brier_score_loss, y_test, test_scores),
        "accuracy_at_threshold": metric_or_none(accuracy_score, y_test, test_pred),
        "precision_at_threshold": metric_or_none(precision_score, y_test, test_pred, zero_division=0),
        "recall_at_threshold": metric_or_none(recall_score, y_test, test_pred, zero_division=0),
        "f1_at_threshold": metric_or_none(f1_score, y_test, test_pred, zero_division=0),
        "notes": [
            "Scores are predicted confidence estimates, not experimentally validated probabilities.",
            "Default training excludes direct curated-evidence fields to reduce label leakage.",
            "Expression features are only informative if the feature table was built with real expression correlations.",
        ],
    }

    rf = pipeline.named_steps["rf"]
    importance_rows = sorted(
        (
            {"feature": feature, "importance": float(importance)}
            for feature, importance in zip(feature_columns, rf.feature_importances_)
        ),
        key=lambda row: row["importance"],
        reverse=True,
    )

    scored_rows = []
    for row, score in zip(rows, all_scores):
        scored = dict(row)
        scored["predicted_confidence"] = f"{score:.6f}"
        scored_rows.append(scored)
    scored_rows.sort(key=lambda row: float(row["predicted_confidence"]), reverse=True)
    for rank, row in enumerate(scored_rows, start=1):
        row["confidence_rank"] = rank

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores_path = args.output_dir / "tf_gene_edge_scores.csv"
    metrics_path = args.output_dir / "edge_confidence_metrics.json"
    importance_path = args.output_dir / "edge_confidence_feature_importance.csv"
    model_path = args.output_dir / "edge_confidence_rf.joblib"

    write_rows(scores_path, ["confidence_rank", "predicted_confidence"] + list(rows[0].keys()), scored_rows)
    write_rows(importance_path, ["feature", "importance"], importance_rows)
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    joblib.dump({
        "pipeline": pipeline,
        "feature_columns": feature_columns,
        "metrics": metrics,
    }, model_path)

    return {
        "scores_path": scores_path,
        "metrics_path": metrics_path,
        "importance_path": importance_path,
        "model_path": model_path,
        "metrics": metrics,
        "top_features": importance_rows[:10],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a random-forest TF-target confidence model.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--min-samples-leaf", type=int, default=2)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument(
        "--include-direct-evidence",
        action="store_true",
        help="Include Evidence/PMID/Source/regulation-role columns. Useful for evidence-completeness scoring, but can leak labels.",
    )
    return parser.parse_args()


def main() -> int:
    result = train_model(parse_args())
    metrics = result["metrics"]
    print(f"Wrote scores: {result['scores_path']}")
    print(f"Wrote metrics: {result['metrics_path']}")
    print(f"Wrote feature importance: {result['importance_path']}")
    print(f"Wrote model: {result['model_path']}")
    print(f"Rows: {metrics['row_count']} | Features: {metrics['feature_count']} | Split: {metrics['split_strategy']}")
    print(f"AUROC: {metrics['roc_auc']} | AUPRC: {metrics['average_precision']} | F1@{metrics['threshold']}: {metrics['f1_at_threshold']}")
    print("Top features:")
    for row in result["top_features"]:
        print(f"  {row['feature']}: {row['importance']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
