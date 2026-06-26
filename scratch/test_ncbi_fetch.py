import urllib.request
import urllib.parse
import json
import xml.etree.ElementTree as ET
import re

def get_promoter_sequence(locus_tag):
    try:
        # 1. Search for the gene ID
        term = f"{locus_tag}[Gene Name] AND 196627[Taxonomy ID]"
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
            "db": "gene",
            "term": term,
            "retmode": "json"
        })
        print(f"Searching for {locus_tag}: {search_url}")
        
        req = urllib.request.Request(search_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            id_list = data.get("esearchresult", {}).get("idlist", [])
            
        if not id_list:
            # Try searching without taxid restriction, just in case
            term = f"{locus_tag}[Locus Tag] AND Corynebacterium glutamicum"
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
                "db": "gene",
                "term": term,
                "retmode": "json"
            })
            with urllib.request.urlopen(urllib.request.Request(search_url, headers={'User-Agent': 'Mozilla/5.0'})) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                id_list = data.get("esearchresult", {}).get("idlist", [])
                
        if not id_list:
            print(f"Gene {locus_tag} not found.")
            return None
            
        gene_id = id_list[0]
        print(f"Found Gene ID: {gene_id}")
        
        # 2. Get gene summary for coordinates
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode({
            "db": "gene",
            "id": gene_id,
            "retmode": "json"
        })
        
        req = urllib.request.Request(summary_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as resp:
            s_data = json.loads(resp.read().decode('utf-8'))
            gene_info = s_data.get("result", {}).get(gene_id, {})
            
        genomic_info = gene_info.get("genomicinfo", [])
        if not genomic_info:
            print(f"No genomic coordinates found for {locus_tag}.")
            return None
            
        g_info = genomic_info[0]
        chr_acc = g_info.get("chraccver")
        chr_start = g_info.get("chrstart")
        chr_stop = g_info.get("chrstop")
        
        if chr_acc is None or chr_start is None or chr_stop is None:
            print(f"Incomplete coordinates for {locus_tag}.")
            return None
            
        print(f"Coordinates: start={chr_start}, stop={chr_stop}, acc={chr_acc}")
        
        # 3. Calculate promoter coordinates (upstream 200bp)
        # Check strand: chrstart > chrstop means negative strand
        is_negative = chr_start > chr_stop
        
        if is_negative:
            # Negative strand: upstream is to the right (higher coordinates)
            prom_start = chr_start + 1
            prom_stop = chr_start + 200
        else:
            # Positive strand: upstream is to the left (lower coordinates)
            prom_start = chr_start - 200
            prom_stop = chr_start - 1
            
        if prom_start < 1:
            prom_start = 1
            
        # 4. Fetch the sequence
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "nuccore",
            "id": chr_acc,
            "seq_start": prom_start,
            "seq_stop": prom_stop,
            "rettype": "fasta",
            "retmode": "text"
        })
        print(f"Fetching sequence: {fetch_url}")
        
        req = urllib.request.Request(fetch_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as resp:
            fasta_data = resp.read().decode('utf-8')
            
        # Extract sequence lines
        lines = fasta_data.strip().splitlines()
        seq_lines = [l.strip() for l in lines if not l.startswith(">")]
        seq = "".join(seq_lines)
        
        # If on the negative strand, we need the reverse complement!
        if is_negative:
            # Reverse complement the sequence
            comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
                    "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}
            rev_comp = "".join(comp.get(base, base) for base in reversed(seq))
            seq = rev_comp
            
        return seq
    except Exception as e:
        print(f"Error fetching promoter for {locus_tag}: {e}")
        return None

# Test with cg0350 (glxR) and cg0986 (amtR)
for locus in ["cg0350", "cg0986"]:
    seq = get_promoter_sequence(locus)
    if seq:
        print(f"Success! {locus} promoter (first 50bp): {seq[:50]}... length={len(seq)}")
    else:
        print(f"Failed to fetch promoter for {locus}")
