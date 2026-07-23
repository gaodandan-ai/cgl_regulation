import csv

import pytest

from scripts.analysis.validate_scientific_outputs import evaluate_gold_standard


def _write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["edge", "score", "label"])
        writer.writeheader()
        writer.writerows(rows)


def test_gold_standard_metrics(tmp_path):
    path = tmp_path / "held_out.csv"
    _write_rows(path, [
        {"edge": "a", "score": 0.9, "label": 1},
        {"edge": "b", "score": 0.8, "label": 1},
        {"edge": "c", "score": 0.2, "label": 0},
        {"edge": "d", "score": 0.1, "label": 0},
    ])
    result = evaluate_gold_standard(path)
    assert result["samples"] == 4
    assert result["metrics"]["auroc"] == 1.0
    assert result["metrics"]["average_precision"] == 1.0
    assert result["metrics"]["brier_score"] == pytest.approx(0.025)
    assert result["metrics"]["expected_calibration_error"] == pytest.approx(0.15)


def test_gold_standard_rejects_single_class(tmp_path):
    path = tmp_path / "invalid.csv"
    _write_rows(path, [
        {"edge": "a", "score": 0.9, "label": 1},
        {"edge": "b", "score": 0.8, "label": 1},
    ])
    with pytest.raises(ValueError, match="both positive and negative"):
        evaluate_gold_standard(path)
