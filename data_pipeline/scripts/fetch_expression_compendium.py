#!/usr/bin/env python3
"""Fetch and summarize the C. glutamicum expression compendium.

The compendium from Kranz, Polen, and Bott contains 927 manually curated
microarray experiments that were combined into 304 microarray sets. This script
downloads the Zenodo files on demand and computes TF-target expression
correlations for the regulatory edges used by this project.

Outputs:
    data/reference/expression_compendium/record_metadata.json
    data/reference/expression_compendium/tf_target_compendium_correlations.csv

The normalized expression workbook is about 91 MB. Use --metadata-only to only
write record metadata and download links.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
OUT_DIR = DATA_DIR / "reference" / "expression_compendium"
RAW_DIR = OUT_DIR / "raw"
ZENODO_RECORD_API = "https://zenodo.org/api/records/6842664"
NORMALIZED_KEY = "Normalized expression matrix_fluorescence intensities.xlsx"
DIFF_EXPR_KEY = "Filtered differential expression results.xlsx"
CORRELATION_OUTPUT = OUT_DIR / "tf_target_compendium_correlations.csv"


def clean(value: Any) -> str:
    return str(value or "").strip()


def normalize_locus(value: Any) -> str:
    raw = clean(value).lower()
    raw = raw.replace("gene:", "")
    match = re.search(r"(?:cg|cgl|ncgl)[_\- ]?(\d{1,4})", raw)
    if match:
        return f"cg{int(match.group(1)):04d}"
    return raw


def load_gene_mapping(path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    alias_to_cg: Dict[str, str] = {}
    cg_to_name: Dict[str, str] = {}
    if not path.exists():
        return alias_to_cg, cg_to_name
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            cg = normalize_locus(row.get("cg_locus"))
            cgl = normalize_locus(row.get("cgl_locus"))
            gene_name = clean(row.get("gene_name"))
            if cg:
                alias_to_cg[cg] = cg
                if gene_name:
                    alias_to_cg[gene_name.lower()] = cg
                    cg_to_name[cg] = gene_name
            if cgl and cg:
                alias_to_cg[cgl] = cg
    return alias_to_cg, cg_to_name


def canonical_gene(value: Any, alias_to_cg: Dict[str, str]) -> str:
    raw = normalize_locus(value)
    if raw in alias_to_cg:
        return alias_to_cg[raw]
    if raw.startswith("cgl"):
        return alias_to_cg.get(raw, raw)
    return raw


def fetch_record() -> Dict[str, Any]:
    request = urllib.request.Request(ZENODO_RECORD_API, headers={"User-Agent": "cgl-regulation/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


def download_file(url: str, path: Path, force: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0 and not force:
        print(f"Using cached file: {path}")
        return
    print(f"Downloading {url} -> {path}")
    request = urllib.request.Request(url, headers={"User-Agent": "cgl-regulation/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)


def read_regulatory_pairs(path: Path, alias_to_cg: Dict[str, str]) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            tf = canonical_gene(row.get("TF_locusTag") or row.get("TF_name"), alias_to_cg)
            target = canonical_gene(row.get("TG_locusTag") or row.get("TG_name"), alias_to_cg)
            if tf and target and tf != target:
                pairs.append((tf, target))
    return sorted(set(pairs))


def first_present(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lower_to_original = {str(col).strip().lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    for col in columns:
        text = str(col).strip().lower()
        if any(candidate.lower() in text for candidate in candidates):
            return col
    return None


def load_expression_matrix(path: Path, alias_to_cg: Dict[str, str]) -> Dict[str, List[float]]:
    try:
        import pandas as pd
    except Exception as exc:
        raise SystemExit(
            "pandas/openpyxl are required to read the compendium workbook. "
            "Install dependencies with: python -m pip install pandas openpyxl"
        ) from exc

    print(f"Reading normalized expression workbook: {path}")
    workbook = pd.ExcelFile(path)
    frames = []
    for sheet in workbook.sheet_names:
        candidate_frames = [
            pd.read_excel(workbook, sheet_name=sheet),
            pd.read_excel(workbook, sheet_name=sheet, header=6),
        ]
        for frame in candidate_frames:
            if frame.empty:
                continue
            gene_columns = [col for col in frame.columns if normalize_locus(col).startswith("cg")]
            if len(gene_columns) >= 100:
                if "Mean/Median" in frame.columns:
                    frame = frame[frame["Mean/Median"].astype(str).str.lower() == "mean"]
                if "Fluorescence" in frame.columns:
                    frame = frame[frame["Fluorescence"].astype(str).str.upper() == "F635"]
                expression: Dict[str, List[float]] = {}
                for col in gene_columns:
                    gene = canonical_gene(col, alias_to_cg)
                    series = pd.to_numeric(frame[col], errors="coerce")
                    values = []
                    for value in series.tolist():
                        try:
                            parsed = float(value)
                        except Exception:
                            parsed = float("nan")
                        if math.isfinite(parsed) and parsed > 0:
                            parsed = math.log2(parsed)
                        else:
                            parsed = float("nan")
                        values.append(parsed)
                    if any(math.isfinite(v) for v in values):
                        expression[gene] = values
                sample_count = len(next(iter(expression.values()))) if expression else 0
                print(
                    f"Loaded expression profiles for {len(expression)} genes across "
                    f"{sample_count} Mean/F635 signal rows from sheet '{sheet}'."
                )
                return expression

        frame = candidate_frames[0]
        if frame.empty:
            continue
        gene_col = first_present(
            list(frame.columns),
            ["cg", "cgl", "gene", "gene_id", "locus", "locus_tag", "identifier", "probe"],
        )
        if not gene_col:
            continue
        frame = frame.rename(columns={gene_col: "gene_id"})
        numeric_cols = [
            col for col in frame.columns
            if col != "gene_id" and getattr(frame[col], "dtype", None) is not None
        ]
        numeric_cols = [col for col in numeric_cols if str(frame[col].dtype).startswith(("int", "float"))]
        if len(numeric_cols) < 3:
            continue
        frames.append(frame[["gene_id", *numeric_cols]])
    if not frames:
        raise SystemExit(f"No usable expression matrix sheet found in {path}")

    combined = frames[0]
    for frame in frames[1:]:
        combined = combined.merge(frame, on="gene_id", how="outer")

    expression: Dict[str, List[float]] = {}
    value_cols = [col for col in combined.columns if col != "gene_id"]
    for _, row in combined.iterrows():
        gene = canonical_gene(row.get("gene_id"), alias_to_cg)
        if not gene:
            continue
        values: List[float] = []
        for col in value_cols:
            value = row.get(col)
            try:
                parsed = float(value)
            except Exception:
                parsed = float("nan")
            values.append(parsed)
        if any(math.isfinite(v) for v in values):
            expression[gene] = values
    print(f"Loaded expression profiles for {len(expression)} genes across {len(value_cols)} columns.")
    return expression


def pearson_with_pvalue(x: Sequence[float], y: Sequence[float]) -> Tuple[Optional[float], Optional[float], int]:
    pairs = [(a, b) for a, b in zip(x, y) if math.isfinite(a) and math.isfinite(b)]
    n = len(pairs)
    if n < 3:
        return None, None, n
    xs = [a for a, _ in pairs]
    ys = [b for _, b in pairs]
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    dx = [a - mean_x for a in xs]
    dy = [b - mean_y for b in ys]
    denom = math.sqrt(sum(a * a for a in dx) * sum(b * b for b in dy))
    if denom == 0:
        return None, None, n
    r = sum(a * b for a, b in zip(dx, dy)) / denom
    r = max(min(r, 1.0), -1.0)
    pvalue = None
    try:
        from scipy import stats
        pvalue = float(stats.pearsonr(xs, ys).pvalue)
    except Exception:
        pvalue = None
    return r, pvalue, n


def write_correlations(
    output: Path,
    pairs: Iterable[Tuple[str, str]],
    expression: Dict[str, List[float]],
    gene_names: Dict[str, str],
) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tf",
        "target",
        "tf_name",
        "target_name",
        "correlation",
        "abs_correlation",
        "pvalue",
        "sample_count",
        "source",
    ]
    count = 0
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for tf, target in pairs:
            tf_values = expression.get(tf)
            target_values = expression.get(target)
            if not tf_values or not target_values:
                continue
            corr, pvalue, n = pearson_with_pvalue(tf_values, target_values)
            if corr is None:
                continue
            writer.writerow({
                "tf": tf,
                "target": target,
                "tf_name": gene_names.get(tf, ""),
                "target_name": gene_names.get(target, ""),
                "correlation": f"{corr:.8g}",
                "abs_correlation": f"{abs(corr):.8g}",
                "pvalue": "" if pvalue is None else f"{pvalue:.8g}",
                "sample_count": n,
                "source": "Kranz2022_expression_compendium_Zenodo6842664",
            })
            count += 1
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch C. glutamicum expression compendium and compute edge correlations.")
    parser.add_argument("--metadata-only", action="store_true", help="Only write Zenodo metadata and download links.")
    parser.add_argument("--download-differential", action="store_true", help="Also cache the smaller filtered differential expression workbook.")
    parser.add_argument("--force-download", action="store_true", help="Re-download cached workbooks.")
    parser.add_argument("--normalized-workbook", type=Path, default=None, help="Use an existing normalized expression workbook.")
    parser.add_argument("--regulations", type=Path, default=DATA_DIR / "reference" / "regulations.csv")
    parser.add_argument("--gene-mapping", type=Path, default=DATA_DIR / "reference" / "gene_mapping.csv")
    parser.add_argument("--output", type=Path, default=CORRELATION_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    record = fetch_record()
    files = {item["key"]: item for item in record.get("files", [])}
    metadata = {
        "title": record.get("metadata", {}).get("title"),
        "doi": record.get("metadata", {}).get("doi"),
        "publication_date": record.get("metadata", {}).get("publication_date"),
        "license": record.get("metadata", {}).get("license", {}).get("id"),
        "files": [
            {"key": item.get("key"), "size": item.get("size"), "url": item.get("links", {}).get("self")}
            for item in record.get("files", [])
        ],
    }
    write_json(OUT_DIR / "record_metadata.json", metadata)
    print(f"Wrote metadata: {OUT_DIR / 'record_metadata.json'}")

    if args.metadata_only:
        return 0

    if args.download_differential and DIFF_EXPR_KEY in files:
        download_file(files[DIFF_EXPR_KEY]["links"]["self"], RAW_DIR / DIFF_EXPR_KEY, args.force_download)

    workbook = args.normalized_workbook
    if workbook is None:
        if NORMALIZED_KEY not in files:
            raise SystemExit(f"Zenodo record does not contain {NORMALIZED_KEY}")
        workbook = RAW_DIR / NORMALIZED_KEY
        download_file(files[NORMALIZED_KEY]["links"]["self"], workbook, args.force_download)

    alias_to_cg, gene_names = load_gene_mapping(args.gene_mapping)
    pairs = read_regulatory_pairs(args.regulations, alias_to_cg)
    expression = load_expression_matrix(workbook, alias_to_cg)
    count = write_correlations(args.output, pairs, expression, gene_names)
    print(f"Wrote {count} TF-target compendium correlations to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
