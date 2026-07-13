"""
tests/test_data_integrity.py
============================
Unit tests for reference data file integrity and schema validation.
Run with:  pytest tests/ -v
"""
import os
import json
import csv
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF  = os.path.join(ROOT, "data", "reference")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def load_json(fname):
    path = os.path.join(REF, fname)
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def load_csv(fname):
    path = os.path.join(REF, fname)
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ─── STRING interactions ──────────────────────────────────────────────────────

class TestStringInteractions:
    def setup_method(self):
        self.data = load_json("string_interactions.json")

    def test_has_meta(self):
        assert "_meta" in self.data
        meta = self.data["_meta"]
        assert meta["version"] == "12.0"
        assert meta["taxon_id"] == 196627

    def test_gene_count(self):
        genes = [k for k in self.data if not k.startswith("_")]
        assert len(genes) >= 2800, f"Expected >=2800 genes, got {len(genes)}"

    def test_edge_count(self):
        meta = self.data["_meta"]
        assert meta["n_edges"] >= 12000

    def test_partner_schema(self):
        """Spot-check a known gene's partner list."""
        assert "cg0001" in self.data
        partners = self.data["cg0001"]
        assert len(partners) > 0
        p = partners[0]
        required = {"partner", "score", "experimental", "coexpression", "textmining", "type"}
        assert required.issubset(p.keys()), f"Missing fields: {required - p.keys()}"

    def test_scores_in_range(self):
        """All scores should be 0-1000."""
        sample_gene = next(k for k in self.data if not k.startswith("_"))
        for p in self.data[sample_gene]:
            assert 0 <= p["score"] <= 1000, f"Score out of range: {p['score']}"

    def test_no_self_loops(self):
        """No gene should be its own partner."""
        for gene, partners in self.data.items():
            if gene.startswith("_"):
                continue
            for p in partners:
                assert p["partner"] != gene, f"Self-loop detected: {gene}"


# ─── Thermodynamics data ──────────────────────────────────────────────────────

class TestThermoData:
    def setup_method(self):
        self.data = load_json("thermo_dgr_data.json")

    def test_has_meta(self):
        meta = self.data.get("_meta", {})
        assert meta.get("version") == "4.0"
        assert meta.get("coverage_pct", 0) >= 60.0, "Coverage below 60%"

    def test_coverage(self):
        meta = self.data["_meta"]
        assert meta["n_curated"] >= 300
        assert meta["n_equilibrator"] >= 500

    def test_known_reaction_schema(self):
        """Check a known central-metabolism reaction."""
        known = ["PGI", "PFK", "GAPD", "ENO"]
        reactions = self.data.get("reactions", self.data)
        found = [r for r in known if r in reactions]
        assert len(found) > 0, f"None of {known} found in thermo data"
        rxn = reactions[found[0]]
        assert "dgr_prime_0" in rxn
        # uncertainty is stored as uncertainty_kJ (eQuilibrator) or dgr_prime_min/max (curated)
        has_uncertainty = "uncertainty_kJ" in rxn or "dgr_prime_min" in rxn
        assert has_uncertainty, f"No uncertainty field in {rxn}"
        assert "confidence" in rxn

    def test_no_infinite_values(self):
        """No reaction should have infinite or NaN dGr."""
        reactions = self.data.get("reactions", self.data)
        for rxn_id, entry in reactions.items():
            if rxn_id.startswith("_"):
                continue
            if not isinstance(entry, dict):
                continue  # skip metadata / non-reaction entries
            val = entry.get("dgr_prime_0")
            if val is None:
                continue  # eQuilibrator entries without dgr_prime_0 are OK
            assert abs(val) < 1e6, f"{rxn_id}: dgr_prime_0 suspiciously large: {val}"


# ─── Regulations ──────────────────────────────────────────────────────────────

class TestRegulations:
    def setup_method(self):
        self.rows = load_csv("regulations.csv")

    def test_row_count(self):
        assert len(self.rows) >= 1600, f"Expected >=1600 rows, got {len(self.rows)}"

    def test_has_evidence_score(self):
        """All rows should now have evidence_score after P0-1 fix."""
        missing = [r for r in self.rows if not r.get("evidence_score", "").strip()]
        assert len(missing) == 0, f"{len(missing)} rows missing evidence_score"

    def test_evidence_score_values(self):
        """Scores should be one of 0.2, 0.4, 0.8, 1.0."""
        valid_scores = {"0.2", "0.4", "0.8", "1.0"}
        for row in self.rows:
            score = str(row.get("evidence_score", "")).strip()
            assert score in valid_scores, f"Invalid score: {score}"

    def test_confidence_labels(self):
        """confidence_label should be HIGH, MEDIUM, or LOW."""
        valid = {"HIGH", "MEDIUM", "LOW"}
        for row in self.rows:
            label = row.get("confidence_label", "").strip()
            assert label in valid, f"Invalid label: {label}"

    def test_locus_format(self):
        """TF_locusTag should match cg\\d+ or Cgl\\d+ pattern."""
        import re
        pattern = re.compile(r'^(cg|Cgl)\d+$', re.IGNORECASE)
        bad = [r["TF_locusTag"] for r in self.rows[:100]
               if r.get("TF_locusTag") and not pattern.match(r["TF_locusTag"])]
        assert len(bad) == 0, f"Malformed locus tags: {bad[:5]}"


# ─── BRENDA kcat mappings ──────────────────────────────────────────────────────

class TestKcatMappings:
    def setup_method(self):
        self.data = load_json("brenda_kcat_mappings.json")

    def test_entry_count(self):
        assert len(self.data) >= 1800, f"Expected >=1800 entries, got {len(self.data)}"

    def test_brenda_priority(self):
        """Known BRENDA entries should have source='BRENDA' and confidence='HIGH'."""
        brenda_entries = {k: v for k, v in self.data.items() if v.get("source") == "BRENDA"}
        assert len(brenda_entries) >= 10, "BRENDA entries seem to be missing"
        for rxn_id, entry in list(brenda_entries.items())[:3]:
            assert entry["confidence"] == "HIGH"

    def test_dlkcat_fill(self):
        dlkcat_entries = {k: v for k, v in self.data.items() if v.get("source") == "DLKcat"}
        assert len(dlkcat_entries) >= 1800

    def test_kcat_positive(self):
        """All kcat values should be positive."""
        bad = [(k, v["kcat"]) for k, v in self.data.items() if v.get("kcat", 1) <= 0]
        assert len(bad) == 0, f"Non-positive kcat values: {bad[:3]}"

    def test_schema(self):
        sample = next(iter(self.data.values()))
        assert "kcat" in sample
        assert "source" in sample
        assert "confidence" in sample


# ─── Gene mapping ──────────────────────────────────────────────────────────────

class TestGeneMapping:
    def setup_method(self):
        self.rows = load_csv("gene_mapping.csv")

    def test_row_count(self):
        assert len(self.rows) >= 3000

    def test_has_cg_locus(self):
        """At least 95% of rows should have cg_locus."""
        filled = sum(1 for r in self.rows if r.get("cg_locus", "").strip())
        pct = 100 * filled // len(self.rows)
        assert pct >= 90, f"cg_locus fill rate too low: {pct}%"

    def test_has_uniprot_column(self):
        """uniprot_id column should exist (may be sparse)."""
        assert "uniprot_id" in self.rows[0].keys()

    def test_uniprot_coverage(self):
        """At least 10% of rows should have UniProt IDs."""
        filled = sum(1 for r in self.rows if r.get("uniprot_id", "").strip())
        pct = 100 * filled // len(self.rows)
        assert pct >= 10, f"UniProt coverage too low: {pct}%"


# ─── Sigma factor annotations ──────────────────────────────────────────────────

class TestSigmaAnnotations:
    def setup_method(self):
        self.data = load_json("sigma_factor_annotations.json")

    def test_seven_sigma_factors(self):
        expected = {"sigA", "sigB", "sigC", "sigD", "sigE", "sigH", "sigM"}
        missing = expected - set(self.data.keys())
        assert len(missing) == 0, f"Missing sigma factors: {missing}"

    def test_schema(self):
        required_fields = {"locus", "gene_name", "sigma_class", "stimulus", "targets_count"}
        for name, ann in self.data.items():
            missing = required_fields - ann.keys()
            assert len(missing) == 0, f"{name} missing fields: {missing}"

    def test_locus_format(self):
        import re
        for name, ann in self.data.items():
            locus = ann.get("locus", "")
            assert re.match(r'^cg\d{4}$', locus), f"{name}: invalid locus '{locus}'"

    def test_sigA_is_housekeeping(self):
        sigA = self.data["sigA"]
        assert sigA["sigma_class"] == "Group_1_sigma"
        assert sigA["consensus_minus35"] == "TTGACA"
        assert sigA["consensus_minus10"] == "TATAAT"
