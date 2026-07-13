"""
backend/sequence_tools.py
=========================
Sequence alignment and homolog search utilities.
"""
import urllib.request
import json

from gene_utils import get_absolute_path

# ── Needleman-Wunsch global alignment ────────────────────────────────────────

def run_needleman_wunsch(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -1):
    n, m   = len(seq1), len(seq2)
    score  = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        score[i][0] = i * gap
    for j in range(m + 1):
        score[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s_match  = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            s_delete = score[i-1][j]   + gap
            s_insert = score[i][j-1]   + gap
            score[i][j] = max(s_match, s_delete, s_insert)

    align1, align2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1; j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] + gap:
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
        else:
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    align1.reverse()
    align2.reverse()
    return "".join(align1), "".join(align2)

# ── Homolog alignment handler ─────────────────────────────────────────────────

_SIMILAR_GROUPS = [
    set("IVLMC"), set("FYW"), set("KR"), set("DE"),
    set("ST"),    set("QN"),  set("AGP"),
]


def handle_homolog_alignment(gene_name: str, accession: str) -> dict:
    if not gene_name or not accession:
        return {"error": "Missing gene_name or accession parameter"}

    # Fetch C. glutamicum sequence
    try:
        req = urllib.request.Request(
            f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            cg_fasta = resp.read().decode("utf-8")
        cg_seq = "".join(cg_fasta.splitlines()[1:])
    except Exception as e:
        return {"error": f"Failed to retrieve sequence for C. glutamicum accession {accession}: {e}"}

    # Find M. tuberculosis homolog
    try:
        for query in (
            f"gene:{gene_name}%20AND%20taxonomy_id:83332",
            f"({gene_name})%20AND%20taxonomy_id:83332",
        ):
            url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&size=1"
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                results = json.loads(resp.read().decode("utf-8")).get("results", [])
            if results:
                break
        if not results:
            return {"error": f"No homolog found in Mycobacterium tuberculosis for gene {gene_name}"}
        homolog_acc  = results[0]["primaryAccession"]
        homolog_org  = results[0]["organism"]["scientificName"]
        homolog_gene = results[0].get("genes", [{}])[0].get("geneName", {}).get("value", gene_name.upper())
    except Exception as e:
        return {"error": f"Failed to search for homolog in M. tuberculosis: {e}"}

    # Fetch M. tuberculosis sequence
    try:
        req = urllib.request.Request(
            f"https://rest.uniprot.org/uniprotkb/{homolog_acc}.fasta",
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            mt_fasta = resp.read().decode("utf-8")
        mt_seq = "".join(mt_fasta.splitlines()[1:])
    except Exception as e:
        return {"error": f"Failed to retrieve sequence for M. tuberculosis accession {homolog_acc}: {e}"}

    # Align
    try:
        a1, a2 = run_needleman_wunsch(cg_seq, mt_seq)
        identity_count = similarity_count = 0
        match_chars = []
        for c1, c2 in zip(a1, a2):
            if c1 == "-" or c2 == "-":
                match_chars.append(" ")
            elif c1 == c2:
                identity_count  += 1
                similarity_count += 1
                match_chars.append("*")
            elif any(c1 in g and c2 in g for g in _SIMILAR_GROUPS):
                similarity_count += 1
                match_chars.append(":")
            else:
                match_chars.append(" ")

        total_len     = len(a1)
        identity_pct  = (identity_count  / total_len * 100) if total_len else 0
        similarity_pct = (similarity_count / total_len * 100) if total_len else 0

        return {
            "cg_accession":          accession,
            "cg_gene_name":          gene_name,
            "homolog_accession":     homolog_acc,
            "homolog_organism":      homolog_org,
            "homolog_gene_name":     homolog_gene,
            "alignment1":            a1,
            "alignment2":            a2,
            "match_string":          "".join(match_chars),
            "identity_percentage":   round(identity_pct,   1),
            "similarity_percentage": round(similarity_pct, 1),
        }
    except Exception as e:
        return {"error": f"Alignment calculation failed: {e}"}
