"""
scripts/fetch_uniprot_ids.py
================================
Fetch UniProt accession IDs for C. glutamicum ATCC 13032 (taxon 196627)
and add uniprot_id + ncbi_protein_id columns to gene_mapping.csv.

Strategy:
  1. Download full proteome from UniProt REST API for taxon 196627
  2. Match by gene_name (exact, case-insensitive)
  3. Match by CGL locus tag (NCgl prefix → Cgl locus)
  4. Write updated gene_mapping.csv

UniProt REST endpoint (no auth required):
  https://rest.uniprot.org/uniprotkb/search?query=organism_id:196627&format=tsv&fields=accession,gene_names,gene_oln,gene_ordered_locus_names,xref_ncbiprotein
"""
import os, csv, urllib.request, urllib.parse, time, re

ROOT       = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MAPPING    = os.path.join(ROOT, "data", "reference", "gene_mapping.csv")
CACHE_FILE = os.path.join(ROOT, "data", "reference", "uniprot_cgl_proteome.tsv")

UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb/search?"
    "query=organism_id%3A196627"
    "&format=tsv"
    "&fields=accession%2Cgene_names%2Cgene_oln%2Cprotein_name%2Creviewed"
    "&size=500"
)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_PARAMS = {
    "query":  "organism_id:196627",
    "format": "tsv",
    "fields": "accession,gene_names,gene_oln,protein_name,reviewed",
    "size":   "500",
}


def download_proteome() -> list:
    """Download all C. glutamicum entries from UniProt (paginated)."""
    if os.path.exists(CACHE_FILE):
        print(f"Using cached proteome: {CACHE_FILE}")
        with open(CACHE_FILE, encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    print("Downloading C. glutamicum proteome from UniProt...")
    all_rows = []
    header   = None

    # Build initial URL
    url = UNIPROT_BASE + "?" + "&".join(f"{k}={urllib.parse.quote(v)}" for k, v in UNIPROT_PARAMS.items())

    page = 0
    while url:
        page += 1
        print(f"  Page {page}: fetching {len(all_rows)} rows so far...")
        req = urllib.request.Request(url, headers={"User-Agent": "cgl_regulation_app/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            # Extract next-page URL from Link header
            link_header = r.headers.get("Link", "")
            next_url = None
            for part in link_header.split(","):
                part = part.strip()
                if 'rel="next"' in part:
                    m = re.search(r'<([^>]+)>', part)
                    if m:
                        next_url = m.group(1)

            content = r.read().decode("utf-8").strip()
            if not content:
                break
            lines = content.split("\n")

            if header is None:
                header = lines[0].split("\t")
                data_lines = lines[1:]
            else:
                data_lines = lines[1:] if lines[0].split("\t")[0] == header[0] else lines

            for line in data_lines:
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) == len(header):
                    all_rows.append(dict(zip(header, fields)))

            url = next_url
        time.sleep(0.3)

    print(f"Downloaded {len(all_rows)} UniProt entries across {page} pages")

    # Cache to TSV
    if all_rows and header:
        with open(CACHE_FILE, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"Cached to {CACHE_FILE}")

    return all_rows


def build_lookup(proteome: list) -> dict:
    """Build multi-key lookup: gene_name/locus -> UniProt accession."""
    lookup = {}  # normalized_key -> {"uniprot_id": str, "ncbi_gene": str}

    for row in proteome:
        acc   = row.get("Entry", "").strip()
        genes = row.get("Gene Names", "").strip()
        oln   = row.get("Gene Names (ordered locus)", "").strip()

        if not acc:
            continue

        entry = {"uniprot_id": acc}

        # Gene names (space-separated, take primary = first)
        for gname in genes.split():
            lookup[gname.lower()] = entry

        # Ordered locus names (e.g. "NCgl1234 Cgl1234")
        for locus in oln.split():
            lookup[locus.lower()] = entry
            # NCgl → cg: NCgl0001 ~ cg0001 numerically
            m = re.match(r'ncgl(\d+)', locus.lower())
            if m:
                num = int(m.group(1))
                cg_tag = f"cg{num:04d}"
                lookup[cg_tag] = entry
            # Cgl → look it up too
            m2 = re.match(r'cgl(\d+)', locus.lower())
            if m2:
                lookup[locus.lower()] = entry

    return lookup


def enrich_mapping(proteome: list):
    """Add uniprot_id column to gene_mapping.csv."""
    lookup = build_lookup(proteome)

    with open(MAPPING, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    fieldnames = list(rows[0].keys())
    if "uniprot_id" not in fieldnames:
        fieldnames.append("uniprot_id")

    n_filled = 0
    for row in rows:
        if row.get("uniprot_id", "").strip():
            continue  # already filled

        matched = None
        # Try cg_locus
        cg = row.get("cg_locus", "").strip().lower()
        if cg and cg in lookup:
            matched = lookup[cg]

        # Try cgl_locus (Cgl0001 format)
        if not matched:
            cgl = row.get("cgl_locus", "").strip().lower()
            if cgl and cgl in lookup:
                matched = lookup[cgl]

        # Try gene_name
        if not matched:
            gname = row.get("gene_name", "").strip().lower()
            if gname and gname in lookup:
                matched = lookup[gname]

        if matched:
            row["uniprot_id"] = matched.get("uniprot_id", "")
            n_filled += 1

    with open(MAPPING, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    print(f"Gene mapping updated: {n_filled}/{total} rows now have UniProt ID ({100*n_filled//total}%)")
    print(f"Written: {MAPPING}")


if __name__ == "__main__":
    proteome = download_proteome()
    print(f"Proteome entries: {len(proteome)}")
    if proteome:
        print("Sample keys:", list(proteome[0].keys()))
        print("Sample row:", proteome[0])
    enrich_mapping(proteome)
