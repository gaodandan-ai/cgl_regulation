#!/usr/bin/env python3
"""Fetch PubMed metadata for PMIDs already cited by local evidence datasets."""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[2]
REF_DIR = ROOT_DIR / "data" / "reference"
OUTPUT_PATH = REF_DIR / "literature" / "pubmed_records.json"


def node_text(node: ET.Element | None) -> str:
    return "" if node is None else "".join(node.itertext()).strip()


def collect_pmids() -> set[str]:
    pmids: set[str] = set()
    for filename in ("regulations.csv", "chipseq_regulations.csv", "regprecise_regulations.csv"):
        with (REF_DIR / filename).open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                pmids.update(re.findall(r"\b\d{6,9}\b", row.get("PMID", "")))
    with (REF_DIR / "tcs_regulations.csv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            pmids.update(re.findall(r"\b\d{6,9}\b", row.get("pmid", "")))
    collectf = REF_DIR / "collectf_tfbs.json"
    if collectf.exists():
        for row in json.loads(collectf.read_text(encoding="utf-8")):
            pmids.update(re.findall(r"\b\d{6,9}\b", str(row.get("pmid", ""))))
    return pmids


def parse_article(article: ET.Element, retrieved_at: str) -> tuple[str, dict] | None:
    pmid = node_text(article.find(".//MedlineCitation/PMID"))
    if not pmid:
        return None
    title = node_text(article.find(".//Article/ArticleTitle"))
    abstract_parts = []
    for element in article.findall(".//Article/Abstract/AbstractText"):
        text = node_text(element)
        label = element.attrib.get("Label")
        if text:
            abstract_parts.append(f"{label}: {text}" if label else text)
    journal = node_text(article.find(".//Article/Journal/Title"))
    year_text = node_text(article.find(".//Article/Journal/JournalIssue/PubDate/Year"))
    if not year_text:
        year_text = node_text(article.find(".//Article/Journal/JournalIssue/PubDate/MedlineDate"))
    year_match = re.search(r"\b(19|20)\d{2}\b", year_text)
    doi = ""
    for article_id in article.findall(".//PubmedData/ArticleIdList/ArticleId"):
        if article_id.attrib.get("IdType") == "doi":
            doi = node_text(article_id)
            break
    return pmid, {
        "title": title,
        "abstract": "\n".join(abstract_parts),
        "journal": journal,
        "year": int(year_match.group(0)) if year_match else None,
        "doi": doi or None,
        "retrieved_at": retrieved_at,
    }


def fetch_pubmed_literature() -> None:
    pmids = sorted(collect_pmids(), key=int)
    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    records: dict[str, dict] = {}
    for start in range(0, len(pmids), 100):
        batch = pmids[start:start + 100]
        query = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "xml"})
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
        request = urllib.request.Request(url, headers={"User-Agent": "cgl-regulation/2.0 literature-indexer"})
        with urllib.request.urlopen(request, timeout=45) as response:
            root = ET.fromstring(response.read())
        for article in root.findall(".//PubmedArticle"):
            parsed = parse_article(article, retrieved_at)
            if parsed:
                records[parsed[0]] = parsed[1]
        if start + 100 < len(pmids):
            time.sleep(0.4)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT_PATH.with_suffix(".json.tmp")
    temporary.write_text(json.dumps({
        "_meta": {"source": "NCBI PubMed E-utilities", "retrieved_at": retrieved_at, "requested": len(pmids), "received": len(records)},
        "records": records,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(OUTPUT_PATH)
    print(f"SUCCESS: cached {len(records)}/{len(pmids)} PubMed records at {OUTPUT_PATH}")


if __name__ == "__main__":
    fetch_pubmed_literature()
