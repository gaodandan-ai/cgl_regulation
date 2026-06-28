import os
import re
import json
import subprocess
import urllib.request
import urllib.parse
import concurrent.futures
from collections import Counter

# 1. Fetch promoter sequence from NCBI (same as working prototype)
def fetch_promoter_single(locus_tag):
    try:
        # Search for the gene ID
        term = f"{locus_tag}[Gene Name] AND 196627[Taxonomy ID]"
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
            "db": "gene",
            "term": term,
            "retmode": "json"
        })
        
        req = urllib.request.Request(search_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            id_list = data.get("esearchresult", {}).get("idlist", [])
            
        if not id_list:
            term = f"{locus_tag}[Locus Tag] AND Corynebacterium glutamicum"
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
                "db": "gene",
                "term": term,
                "retmode": "json"
            })
            with urllib.request.urlopen(urllib.request.Request(search_url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=5) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                id_list = data.get("esearchresult", {}).get("idlist", [])
                
        if not id_list:
            return None
            
        gene_id = id_list[0]
        
        # Get coordinates
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode({
            "db": "gene",
            "id": gene_id,
            "retmode": "json"
        })
        
        req = urllib.request.Request(summary_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            s_data = json.loads(resp.read().decode('utf-8'))
            gene_info = s_data.get("result", {}).get(gene_id, {})
            
        genomic_info = gene_info.get("genomicinfo", [])
        if not genomic_info:
            return None
            
        g_info = genomic_info[0]
        chr_acc = g_info.get("chraccver")
        chr_start = g_info.get("chrstart")
        chr_stop = g_info.get("chrstop")
        
        if chr_acc is None or chr_start is None or chr_stop is None:
            return None
            
        is_negative = chr_start > chr_stop
        if is_negative:
            prom_start = chr_start + 1
            prom_stop = chr_start + 200
        else:
            prom_start = chr_start - 200
            prom_stop = chr_start - 1
            
        if prom_start < 1:
            prom_start = 1
            
        # Fetch sequence
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "nuccore",
            "id": chr_acc,
            "seq_start": prom_start,
            "seq_stop": prom_stop,
            "rettype": "fasta",
            "retmode": "text"
        })
        
        req = urllib.request.Request(fetch_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            fasta_data = resp.read().decode('utf-8')
            
        lines = fasta_data.strip().splitlines()
        seq_lines = [l.strip() for l in lines if not l.startswith(">")]
        seq = "".join(seq_lines).upper()
        
        if is_negative:
            comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
            seq = "".join(comp.get(base, base) for base in reversed(seq))
            
        return seq
    except Exception as e:
        print(f"Error fetching {locus_tag}: {e}")
        return None

def fetch_promoters_parallel(genes):
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_gene = {executor.submit(fetch_promoter_single, g): g for g in genes}
        for future in concurrent.futures.as_completed(future_to_gene):
            gene = future_to_gene[future]
            try:
                seq = future.result()
                if seq:
                    results[gene] = seq
            except Exception as e:
                print(f"Gene {gene} generated an exception: {e}")
    return results

# 2. Fallback Motif Finder: Find consensus k-mer and construct PWM
def find_motif_fallback(sequences, k=10):
    if not sequences:
        return None
    
    # Extract all k-mers from all sequences
    kmers = []
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if "N" not in kmer:
                kmers.append(kmer)
                
    if not kmers:
        return None
        
    # Count k-mers
    kmer_counts = Counter(kmers)
    
    # We want to find a k-mer that is representative, allowing 1 mismatch
    def get_hamming_distance(s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
        
    # To be fast, let's look at the top 100 most frequent k-mers and score them
    # by how many sequences contain them (allowing 1 mismatch)
    top_candidates = [item[0] for item in kmer_counts.most_common(100)]
    best_candidate = None
    best_score = -1
    best_matches = []
    
    for candidate in top_candidates:
        matches = []
        for seq in sequences:
            # Find best match in this sequence
            best_seq_match = None
            min_dist = 999
            for i in range(len(seq) - k + 1):
                sub = seq[i:i+k]
                if "N" in sub:
                    continue
                dist = get_hamming_distance(candidate, sub)
                if dist < min_dist:
                    min_dist = dist
                    best_seq_match = sub
            # If distance <= 2, we consider it a match
            if min_dist <= 2:
                matches.append(best_seq_match)
        
        # Score is number of sequences matched
        score = len(matches)
        if score > best_score:
            best_score = score
            best_candidate = candidate
            best_matches = matches
            
    if not best_candidate or not best_matches:
        return None
        
    # Build Position Weight Matrix (PWM) from matching sites
    pwm = []
    for col in range(k):
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        for match in best_matches:
            char = match[col]
            if char in counts:
                counts[char] += 1
        total = sum(counts.values()) or 1
        # Add pseudocount of 0.1 to avoid 0 probabilities
        pwm.append({
            "A": (counts["A"] + 0.1) / (total + 0.4),
            "C": (counts["C"] + 0.1) / (total + 0.4),
            "G": (counts["G"] + 0.1) / (total + 0.4),
            "T": (counts["T"] + 0.1) / (total + 0.4),
        })
        
    # Generate consensus sequence
    bases = ["A", "C", "G", "T"]
    consensus = "".join(max(bases, key=lambda b: pos[b]) for pos in pwm)
    
    return {
        "consensus": consensus,
        "pwm": pwm,
        "nsites": len(best_matches)
    }

# Test fallback motif finder on simulated sequences
sim_seqs = [
    "AGTCATGTGACTGTTCACACAGTC",
    "TTTGTGACTGTTCACACCCT",
    "GCTGTGACATATCACAGGCG",
    "AATGTGACATCTCACATTTT"
]
# They all contain a variation of "TGTGACNNNTCACA" (GlxR motif)
res = find_motif_fallback(sim_seqs, k=14)
print("Fallback consensus:", res["consensus"])
print("First position PWM:", res["pwm"][0])
