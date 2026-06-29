import math

def get_reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def calculate_window_affinity(pwm: list, window: str, temp_celsius: float = 30.0) -> float:
    # Gas constant R in kcal/(mol * K)
    R = 1.987e-3
    T = temp_celsius + 273.15
    RT = R * T
    
    # Reference binding energy for consensus sequence
    delta_G_ref = -10.0  # kcal/mol (standard strong affinity)
    
    # Compute consensus probabilities
    consensus_probs = []
    for pos_probs in pwm:
        max_base = max(pos_probs, key=pos_probs.get)
        consensus_probs.append(max(1e-4, pos_probs[max_base]))
        
    # Compute energy penalty for mismatch from consensus
    energy_penalty = 0.0
    for j, base in enumerate(window):
        if j >= len(pwm):
            break
        pos_probs = pwm[j]
        base_prob = max(1e-4, pos_probs.get(base, 0.25))
        cons_prob = consensus_probs[j]
        
        # Energy deviation penalty: delta_delta_G = -RT * ln(P(base)/P(consensus))
        delta_delta_G = -RT * math.log(base_prob / cons_prob)
        energy_penalty += delta_delta_G
        
    delta_G_bind = delta_G_ref + energy_penalty
    return delta_G_bind

def scan_sequence_for_affinity(pwm: list, sequence: str, temp_celsius: float = 30.0) -> dict:
    L = len(pwm)
    N = len(sequence)
    if N < L:
        return {"error": f"Sequence length ({N}) is shorter than PWM motif width ({L})."}
        
    R = 1.987e-3
    T = temp_celsius + 273.15
    RT = R * T
    
    best_delta_G = 999.0
    best_pos = -1
    best_strand = "+"
    best_seq = ""
    
    sequence_upper = sequence.upper()
    profile = []
    
    for i in range(N - L + 1):
        window_fwd = sequence_upper[i : i + L]
        dG_fwd = calculate_window_affinity(pwm, window_fwd, temp_celsius)
        
        window_rev = get_reverse_complement(window_fwd)
        dG_rev = calculate_window_affinity(pwm, window_rev, temp_celsius)
        
        if dG_fwd <= dG_rev:
            pos_dG = dG_fwd
            pos_strand = "+"
            pos_seq = window_fwd
        else:
            pos_dG = dG_rev
            pos_strand = "-"
            pos_seq = window_rev
            
        pos_kd = 1e9 * math.exp(pos_dG / RT) if pos_dG / RT < 50 else 1e9
        
        profile.append({
            "position": i + 1,
            "strand": pos_strand,
            "matched_sequence": pos_seq,
            "delta_G": pos_dG,
            "Kd": pos_kd
        })
        
        if pos_dG < best_delta_G:
            best_delta_G = pos_dG
            best_pos = i + 1
            best_strand = pos_strand
            best_seq = pos_seq
            
    # Compute Kd in nM
    exponent = best_delta_G / RT
    if exponent > 50:
        kd_nm = 1e9
    elif exponent < -50:
        kd_nm = 1e-9
    else:
        kd_nm = 1e9 * math.exp(exponent)
        
    return {
        "position": best_pos,
        "strand": best_strand,
        "matched_sequence": best_seq,
        "delta_G": best_delta_G,
        "Kd": kd_nm,
        "profile": profile
    }
