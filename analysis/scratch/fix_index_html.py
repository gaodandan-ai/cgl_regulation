# fix_index_html.py
with open('web/index.html', 'r', encoding='utf-8') as f:
    content = f.read()

# 1. Revert left sidebar changes (remove the protein-analysis-section block)
# Let's find: <!-- 蛋白结构与结合分析 --> or similar
# Let's locate starting from <!-- 蛋白结构与结合分析 --> down to the next action-section comment
revert_start = content.find("<!-- 蛋白结构与结合分析 -->")
if revert_start == -1:
    # try matching with whatever is on line 310
    revert_start = content.find("<!-- 蛋白结构与结合分析") # fallback substring
if revert_start == -1:
    revert_start = content.find("<!-- 蛋白结构与结合")
if revert_start == -1:
    # Let's find based on class name
    revert_start = content.find('<div class="sidebar-section protein-analysis-section">')
    # if so, back track to comments
    if revert_start != -1:
        revert_start = content.rfind("<!--", 0, revert_start)

revert_end = content.find('<div class="sidebar-section action-section">', revert_start)

if revert_start != -1 and revert_end != -1:
    print(f"Reverting left sidebar block from {revert_start} to {revert_end}...")
    content = content[:revert_start] + content[revert_end:]
    print("Revert Succeeded!")
else:
    print("Could not find left sidebar block to revert!")

# 2. Insert the two new sections into the right details panel
# Locate the end of the first detail-section (基本信息) which ends with </table>\\n                    </div>
# Let's search for "info-pathway-container" or "info-table" or the closing div of "基本信息" detail-section
basic_info_close = content.find("</table>\n                    </div>", content.find("info-locus"))
if basic_info_close == -1:
    basic_info_close = content.find("</table>\r\n                    </div>", content.find("info-locus"))
if basic_info_close == -1:
    basic_info_close = content.find("</table>", content.find("info-locus"))
    if basic_info_close != -1:
        # find closing div after table
        basic_info_close = content.find("</div>", basic_info_close) + len("</div>")
else:
    basic_info_close += len("</table>\n                    </div>") # move index past the end of the closing div

print("Insert index basic_info_close:", basic_info_close)

if basic_info_close != -1:
    new_sections = """

                    <!-- 蛋白结构域 & 结合基序预测 -->
                    <div class="detail-section protein-analysis-section" id="detail-protein-domain-section" style="display: none;">
                        <h3><i class="fa-solid fa-shapes" style="color: #7c3aed;"></i> 蛋白结构域 & 结合基序预测</h3>
                        <div class="protein-predict-container" style="display: grid; grid-template-columns: 90px 1fr; gap: 10px; align-items: start; margin-bottom: 8px;">
                            <!-- 左列：蛋白质三维卡通结构图 -->
                            <div class="protein-structure-col" style="text-align: center;">
                                <img src="lib/protein_structure_mock.png" alt="Protein Structure" style="width: 80px; height: 80px; object-fit: contain; border-radius: 6px; border: 1px solid var(--border-color); background: #ffffff; padding: 4px; box-shadow: 0 2px 8px rgba(0,0,0,0.04);">
                                <div style="font-size: 9px; color: var(--text-muted); margin-top: 4px; line-height: 1.2; font-weight: 500;">三维结构模型</div>
                            </div>
                            <!-- 右列：WebLogo + 热图 -->
                            <div class="protein-motif-col" style="display: flex; flex-direction: column; gap: 4px;">
                                <div style="font-size: 10px; color: var(--text-secondary); font-weight: 600; display: flex; align-items: center; gap: 4px;"><i class="fa-solid fa-align-left" style="font-size: 8px; color: #7c3aed;"></i> 结合基序 WebLogo:</div>
                                <canvas id="right-motif-logo-canvas" width="220" height="50" style="background:#ffffff; border:1px solid var(--border-color); border-radius:4px; display:block; width:100%; height:45px;"></canvas>
                                
                                <div style="font-size: 10px; color: var(--text-secondary); font-weight: 600; margin-top: 4px; display: flex; align-items: center; gap: 4px;"><i class="fa-solid fa-border-all" style="font-size: 8px; color: #7c3aed;"></i> 位置权重矩阵 (PWM) 热图:</div>
                                <canvas id="right-motif-heatmap-canvas" width="220" height="55" style="background:#ffffff; border:1px solid var(--border-color); border-radius:4px; display:block; width:100%; height:50px;"></canvas>
                            </div>
                        </div>
                        <div id="right-protein-domain-result" style="margin-top: 8px; font-size: 11px; color: var(--text-secondary); line-height: 1.45; max-height: 150px; overflow-y: auto; background: var(--bg-card); padding: 8px 10px; border-radius: 6px; border: 1px solid var(--border-color);"></div>
                    </div>

                    <!-- 蛋白结合位点 & 占位分析 -->
                    <div class="detail-section binding-analysis-section" id="detail-binding-site-section" style="display: none;">
                        <h3><i class="fa-solid fa-chart-line" style="color: #dc2626;"></i> 蛋白结合位点 & 占位分析</h3>
                        <p style="font-size: 11px; color: var(--text-secondary); margin-bottom: 8px; line-height: 1.45;">
                            ChIP-seq 结合丰度 (染色体区域结合强度)。波峰越高，结合富集强度越强。
                        </p>
                        <!-- ChIP-seq 占位彩虹图 -->
                        <canvas id="right-chipseq-peak-canvas" width="310" height="90" style="background:#ffffff; border:1px solid var(--border-color); border-radius:6px; margin-bottom:8px; display:block; width:100%; height:85px;"></canvas>
                        
                        <!-- 占位条件切换 -->
                        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:8px; background:rgba(99, 102, 241, 0.04); padding:5px 10px; border-radius:6px; border:1px solid rgba(99, 102, 241, 0.1);">
                            <span style="font-size:9.5px; color:var(--text-secondary); font-weight:600;"><i class="fa-solid fa-sliders"></i> 环境响应模拟:</span>
                            <div style="display:flex; gap:4px;">
                                <button id="btn-right-cond-ctrl" class="secondary-btn active" style="font-size:8.5px; padding:2px 6px; height:18px; border-radius:3px; width:auto; border:1px solid var(--color-primary-accent); background-color:rgba(30, 58, 138, 0.08); color:var(--color-primary-accent); font-weight:600; cursor:pointer;">Control</button>
                                <button id="btn-right-cond-stress" class="secondary-btn" style="font-size:8.5px; padding:2px 6px; height:18px; border-radius:3px; width:auto; border:1px solid var(--border-color); background-color:#ffffff; color:var(--text-secondary); cursor:pointer;">Stress</button>
                            </div>
                        </div>

                        <!-- 结合位点表格 -->
                        <div style="margin-top: 8px; overflow-x: auto; border: 1px solid var(--border-color); border-radius: 6px;">
                            <table class="relations-table" id="right-binding-sites-table" style="width: 100%; font-size: 10.5px; border-collapse: collapse;">
                                <thead>
                                    <tr style="background-color: var(--bg-card); border-bottom: 1px solid var(--border-color);">
                                        <th style="padding: 6px 8px; text-align: left; font-weight: 600; color: var(--text-secondary);">DNA 结合位点</th>
                                        <th style="padding: 6px 8px; text-align: left; font-weight: 600; color: var(--text-secondary);">基因组位置</th>
                                        <th style="padding: 6px 8px; text-align: right; font-weight: 600; color: var(--text-secondary);">ChIP-seq occupancy</th>
                                    </tr>
                                </thead>
                                <tbody style="font-family: monospace;">
                                    <!-- 动态渲染 -->
                                </tbody>
                            </table>
                        </div>
                    </div>"""
                    
    content = content[:basic_info_close] + new_sections + content[basic_info_close:]
    with open('web/index.html', 'w', encoding='utf-8') as f:
        f.write(content)
    print("Injection Succeeded!")
else:
    print("Could not find insert index!")
