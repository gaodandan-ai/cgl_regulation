# fix_app_js.py
with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

# Let's target the exact text from the previous output:
target = """            let htmlContent = '';

            if (result.total_sites !== undefined) {

                htmlContent += `<div style="font-size: 11px; font-weight: 600; color: var(--color-primary-accent); margin-bottom: 8px;"><i class="fa-solid fa-list-check"></i> ػǼǽλ: ${result.total_sites}</div>`;

            }

            htmlContent += parseMarkdownToHtml(result.summary || '޷');


            resultCard.innerHTML = htmlContent;"""

# Let's inspect if target is in content
# Since characters might have encoding issues, let's look for a substring first
sub_target = "resultCard.classList.remove('loading');"
idx = content.find(sub_target)
print("Index of sub_target:", idx)

# We will replace the block from resultCard.classList.remove('loading'); down to resultCard.innerHTML = htmlContent;
# Let's find the start and end of this block:
block_start = content.find("resultCard.classList.remove('loading');", idx)
block_end = content.find("resultCard.innerHTML = htmlContent;", block_start) + len("resultCard.innerHTML = htmlContent;")

print("Start:", block_start)
print("End:", block_end)

if block_start != -1 and block_end != -1:
    replacement = """resultCard.classList.remove('loading');

            // Look up local sequences for motif prediction
            const tfLower = query.toLowerCase();
            const localSeqs = [];
            regulations.forEach(row => {
                if ((row.TF_locusTag && row.TF_locusTag.toLowerCase() === tfLower) || 
                    (row.TF_name && row.TF_name.toLowerCase() === tfLower)) {
                    const site = row.Binding_site;
                    if (site && site.trim() && site.trim() !== 'nan') {
                        site.split(';').forEach(s => {
                            if (s.trim()) localSeqs.push(s.trim());
                        });
                    }
                }
            });

            let htmlContent = '';
            if (result.total_sites !== undefined) {
                htmlContent += `<div style="font-size: 11.5px; font-weight: 600; color: var(--color-primary-accent); margin-bottom: 8px;"><i class="fa-solid fa-list-check"></i> 本地基因组已知结合靶标: ${result.total_sites} 个</div>`;
            }

            // Add visual canvases for motif logo and ChIP-seq signal track
            htmlContent += `
                <div class="visual-analysis-panel" style="margin-top: 12px; margin-bottom: 12px;">
                    <div style="font-size: 11px; font-weight: 600; color: var(--text-secondary); margin-bottom: 6px;"><i class="fa-solid fa-shapes"></i> 结合基序预测 (Binding Motif Logo):</div>
                    <canvas id="motif-logo-canvas" width="280" height="70" style="background:#ffffff; border:1px solid var(--border-color); border-radius:6px; margin-bottom:10px; width:100%; display:block;"></canvas>

                    <div style="font-size: 11px; font-weight: 600; color: var(--text-secondary); margin-bottom: 6px;"><i class="fa-solid fa-chart-area"></i> ChIP-seq 占位信号丰度 (Signal Intensity):</div>
                    <canvas id="chipseq-peak-canvas" width="280" height="110" style="background:#ffffff; border:1px solid var(--border-color); border-radius:6px; margin-bottom:8px; width:100%; display:block;"></canvas>
                    
                    <div style="display:flex; justify-content:space-between; align-items:center; background:rgba(99, 102, 241, 0.04); padding:6px 10px; border-radius:6px; border:1px solid rgba(99, 102, 241, 0.1);">
                        <span style="font-size:10px; color:var(--text-secondary); font-weight:600;"><i class="fa-solid fa-sliders"></i> 环境响应模拟:</span>
                        <div style="display:flex; gap:6px;">
                            <button id="btn-cond-ctrl" class="secondary-btn active" style="font-size:9px; padding:2px 8px; height:20px; border-radius:4px; width:auto; border:1px solid var(--color-primary-accent); background-color:rgba(30, 58, 138, 0.08); color:var(--color-primary-accent); font-weight:600;">Control</button>
                            <button id="btn-cond-stress" class="secondary-btn" style="font-size:9px; padding:2px 8px; height:20px; border-radius:4px; width:auto; border:1px solid var(--border-color); background-color:#ffffff; color:var(--text-secondary);">Stress</button>
                        </div>
                    </div>
                </div>
            `;

            htmlContent += parseMarkdownToHtml(result.summary || '无分析数据');
            resultCard.innerHTML = htmlContent;

            // Render Motif Logo
            const motifCanvas = document.getElementById('motif-logo-canvas');
            if (motifCanvas) {
                renderMotifLogo(motifCanvas, localSeqs);
            }

            // Render ChIP-seq peak with animations/toggles
            const chipseqCanvas = document.getElementById('chipseq-peak-canvas');
            const btnCtrl = document.getElementById('btn-cond-ctrl');
            const btnStress = document.getElementById('btn-cond-stress');

            if (chipseqCanvas) {
                let currentScale = 0.75;
                let currentCond = 'Control';
                renderChipSeqPeak(chipseqCanvas, currentScale, currentCond);

                if (btnCtrl && btnStress) {
                    btnCtrl.addEventListener('click', () => {
                        btnCtrl.classList.add('active');
                        btnCtrl.style.borderColor = 'var(--color-primary-accent)';
                        btnCtrl.style.backgroundColor = 'rgba(30, 58, 138, 0.08)';
                        btnCtrl.style.color = 'var(--color-primary-accent)';
                        btnCtrl.style.fontWeight = '600';

                        btnStress.classList.remove('active');
                        btnStress.style.borderColor = 'var(--border-color)';
                        btnStress.style.backgroundColor = '#ffffff';
                        btnStress.style.color = 'var(--text-secondary)';
                        btnStress.style.fontWeight = '500';

                        currentScale = 0.75;
                        currentCond = 'Control';
                        renderChipSeqPeak(chipseqCanvas, currentScale, currentCond);
                    });

                    btnStress.addEventListener('click', () => {
                        btnStress.classList.add('active');
                        btnStress.style.borderColor = '#ef4444';
                        btnStress.style.backgroundColor = 'rgba(239, 68, 68, 0.08)';
                        btnStress.style.color = '#ef4444';
                        btnStress.style.fontWeight = '600';

                        btnCtrl.classList.remove('active');
                        btnCtrl.style.borderColor = 'var(--border-color)';
                        btnCtrl.style.backgroundColor = '#ffffff';
                        btnCtrl.style.color = 'var(--text-secondary)';
                        btnCtrl.style.fontWeight = '500';

                        currentScale = 0.25;
                        currentCond = 'Stress';
                        renderChipSeqPeak(chipseqCanvas, currentScale, currentCond);
                    });
                }
            }"""
            
    new_content = content[:block_start] + replacement + content[block_end:]
    with open('web/app.js', 'w', encoding='utf-8') as f:
        f.write(new_content)
    print("Replacement Succeeded!")
else:
    print("Block not found!")
