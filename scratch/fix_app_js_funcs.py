# fix_app_js_funcs.py
with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

# Let's append the functions at the end of the file
additional_code = """

// ==========================================================================
// 11. Visual Chart Drawing Helpers (Motif Logo & ChIP-seq Peak)
// ==========================================================================

function renderMotifLogo(canvas, sequences) {
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    ctx.clearRect(0, 0, width, height);

    let pwm = [];
    let motifLen = 12;
    
    if (sequences && sequences.length > 0) {
        const cleanSeqs = [];
        let maxLen = 0;
        sequences.forEach(s => {
            const seq = s.trim().toUpperCase().replace(/[^ACGT]/g, '');
            if (seq.length > 0) {
                cleanSeqs.push(seq);
                if (seq.length > maxLen) maxLen = seq.length;
            }
        });
        
        if (cleanSeqs.length > 0) {
            motifLen = Math.min(16, maxLen);
            for (let i = 0; i < motifLen; i++) {
                const counts = { A: 0, C: 0, G: 0, T: 0 };
                cleanSeqs.forEach(seq => {
                    if (i < seq.length) {
                        const char = seq[i];
                        if (counts[char] !== undefined) counts[char]++;
                    }
                });
                const total = Object.values(counts).reduce((a, b) => a + b, 0) || 1;
                pwm.push({
                    A: counts.A / total,
                    C: counts.C / total,
                    G: counts.G / total,
                    T: counts.T / total
                });
            }
        }
    }
    
    if (pwm.length === 0) {
        const mockSeq = "TGTGACGTGTCT";
        motifLen = mockSeq.length;
        for (let i = 0; i < motifLen; i++) {
            const char = mockSeq[i];
            const freqs = { A: 0.05, C: 0.05, G: 0.05, T: 0.05 };
            freqs[char] = 0.85;
            pwm.push(freqs);
        }
    }

    const colWidth = width / motifLen;
    ctx.textBaseline = 'bottom';
    
    for (let pos = 0; pos < motifLen; pos++) {
        const freqs = pwm[pos];
        const sortedBases = Object.entries(freqs).sort((a, b) => a[1] - b[1]);
        
        let currentY = height - 12;
        const availableHeight = height - 18;
        
        sortedBases.forEach(([base, freq]) => {
            if (freq < 0.02) return;
            const letterHeight = freq * availableHeight;
            
            ctx.save();
            ctx.font = `bold ${Math.max(10, Math.round(letterHeight * 1.3))}px 'Outfit', sans-serif`;
            
            if (base === 'A') ctx.fillStyle = '#2e7d32';
            else if (base === 'C') ctx.fillStyle = '#1976d2';
            else if (base === 'G') ctx.fillStyle = '#e65100';
            else if (base === 'T') ctx.fillStyle = '#d32f2f';
            
            ctx.translate(pos * colWidth + colWidth / 2, currentY);
            ctx.scale(colWidth / 18, letterHeight / 20);
            ctx.fillText(base, -8, 0);
            ctx.restore();
            
            currentY -= letterHeight;
        });
    }
    
    ctx.font = '8px monospace';
    ctx.fillStyle = '#64748b';
    ctx.textAlign = 'center';
    for (let pos = 0; pos < motifLen; pos++) {
        ctx.fillText(pos + 1, pos * colWidth + colWidth / 2, height - 1);
    }
}

function renderChipSeqPeak(canvas, occupancyScale, conditionName) {
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    const width = canvas.width;
    const height = canvas.height;
    
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = '#f8fafc';
    ctx.fillRect(0, 0, width, height);
    
    ctx.strokeStyle = '#e2e8f0';
    ctx.lineWidth = 1;
    
    ctx.beginPath();
    ctx.moveTo(0, height - 20); ctx.lineTo(width, height - 20);
    ctx.moveTo(width / 2, 0); ctx.lineTo(width / 2, height - 20);
    ctx.stroke();
    
    const peakHeight = (height - 35) * occupancyScale;
    const peakCenter = width * 0.4;
    
    const grad = ctx.createLinearGradient(0, height - 20 - peakHeight, 0, height - 20);
    if (conditionName === 'Stress') {
        grad.addColorStop(0, 'rgba(239, 68, 68, 0.7)');
        grad.addColorStop(1, 'rgba(239, 68, 68, 0.05)');
        ctx.strokeStyle = '#ef4444';
    } else {
        grad.addColorStop(0, 'rgba(30, 58, 138, 0.6)');
        grad.addColorStop(1, 'rgba(30, 58, 138, 0.05)');
        ctx.strokeStyle = '#1e3a8a';
    }
    
    ctx.fillStyle = grad;
    ctx.lineWidth = 2;
    
    ctx.beginPath();
    ctx.moveTo(0, height - 20);
    
    ctx.bezierCurveTo(
        peakCenter - 60, height - 20,
        peakCenter - 30, height - 20 - peakHeight,
        peakCenter, height - 20 - peakHeight
    );
    ctx.bezierCurveTo(
        peakCenter + 30, height - 20 - peakHeight,
        peakCenter + 60, height - 20,
        width, height - 20
    );
    ctx.lineTo(width, height - 20);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
    
    ctx.fillStyle = '#475569';
    ctx.font = '8px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText(`条件: ${conditionName === 'Stress' ? '逆境胁迫 (Stress)' : '对照组 (Control)'}`, 8, 14);
    
    ctx.textAlign = 'right';
    ctx.fillText(`信号丰度 (Occupancy): ${Math.round(occupancyScale * 100)}%`, width - 8, 14);
    
    ctx.strokeStyle = '#64748b';
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    ctx.moveTo(width / 2, 25);
    ctx.lineTo(width / 2, height - 20);
    ctx.stroke();
    ctx.setLineDash([]);
    
    ctx.fillStyle = '#64748b';
    ctx.font = '8px monospace';
    ctx.textAlign = 'center';
    ctx.fillText('TSS (0)', width / 2, 22);
    
    ctx.fillText('-200bp', 24, height - 8);
    ctx.fillText('-35bp', peakCenter, height - 8);
    ctx.fillText('+50bp', width - 24, height - 8);
}
"""

with open('web/app.js', 'a', encoding='utf-8') as f:
    f.write(additional_code)

print("Visual functions successfully appended!")
