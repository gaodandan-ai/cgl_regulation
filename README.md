# Cgl Regulation Explorer (v0.5.0)

Bilingual Version: [English](#english-version) | [中文版本](#中文版本)

---

## English Version



An interactive, reproducible regulatory-metabolic analysis platform for **Corynebacterium glutamicum** DSM 20300 / ATCC 13032.

**Live Production URL**: [https://cgl-regulation.vercel.app/](https://cgl-regulation.vercel.app/)

Cgl Regulation Explorer connects transcription factor (TF) and sRNA regulatory network evidence with genome-scale metabolic reaction equations and protein constraints to support hypothesis generation and strain engineering target prioritization.

---

### 1. Platform Workflow Pipeline

The platform processes data across 8 distinct analytical phases:

1. **Regulatory Evidence**: Integrates curated databases, binding motifs, ChIP-seq records, and sRNA-mRNA prediction profiles.
2. **RF Confidence Scoring**: Computes prioritization scores using a Random Forest machine learning model trained on evidence feature matrices.
3. **Metabolic Model Mapping**: Maps target genes to reaction equations and pathways in the *iCW773* genome-scale model.
4. **ecCGL1 Enzyme Annotations**: Enriches reactions with protein context, including kcat values, molecular weight (MW), EC numbers, and parent reactions.
5. **Pathway-centered Analysis**: Traces likely upstream regulators controlling active pathway modules.
6. **Prioritized Engineering Targets**: Ranks candidate regulators by metabolic reaction coverage and confidence.
7. **Systemic Topology Roles**: Annotates global regulators (Global Hubs, Modular Regulators, and Pathway Genes) based on **Abasy Atlas** topological properties.
8. **Case Study Reports**: Generates reproducible case reports and downloadable JSON summaries (e.g. Glutamate-associated regulation, TCA cycle upstream regulators, and Amino Acid biosyntheses).

---

### 2. Main Features

- **Gene / TF Explorer**: Query target networks, visual evidence, operon context, and metabolic mappings.
- **Thermodynamic Promoter Scanner**: Interactive sliding-window DNA-binding energy scanning using position weight matrices (PWMs) and Chart.js line visualizations.
- **Pathway View**: Analyze pathway-specific upstream TFs and reaction lists, accompanied by interactive Cytoscape.js compound-reaction flow diagrams.
- **Active Simulation & Biophysical Affinities**:
  - Runs dynamic regulatory flux balance analysis and enzyme-constrained FBA (ecFBA) locally or via serverless fallbacks.
  - **Enzyme Bottlenecks Analyzer**: Pinpoints metabolic bottlenecks, shadow prices, and capacity constraints under varying temperatures (scaled via Arrhenius equation).
- **Network Topology Section**: Calculates graph metrics (centrality, degrees) and color-codes nodes dynamically by pleiotropic risk (red for global hubs, blue for modular regulators, green for pathway genes) derived from Abasy Atlas.
- **iModulon Explorer**: Integrates independent components (iModulons) to cluster genes by co-expression under varying growth conditions.
- **Engineering Targets Page**: Rank TFs globally or by specific pathway keyword filters using candidate scores and Global Metabolic Impact rankings.
- **Data & Model Quality Dashboard**: Audits platform-wide node/edge counts, edge confidence ranges, and model mapping stats.
- **RAG-based AI Summary**: Integrates local literature databases with LLMs (Gemini, OpenAI, DeepSeek, etc.) to generate automated functional summaries for queried genes.
- **Data Exports**: Download network maps, table records, and audit metrics as PNG, CSV, or JSON.

---

### 3. Data Sources & References

Please cite the original databases and publications when utilizing derived platform results:

| Source Database / Model          | Description                                                                                        | Citations & Reference Links                                                                                                                 |
| :------------------------------- | :------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------ |
| **CoryneRegNet 7**         | Reference database for corynebacterial gene regulatory networks (TFs, binding sites, and targets). | Parise, M.T.D. et al.*Scientific Data* 7, 142 (2020). [DOI: 10.1038/s41597-020-0481-0](https://doi.org/10.1038/s41597-020-0481-0)          |
| **Abasy Atlas**            | Systemic topological roles mapping (Global Hubs, Modular Regulators, Pathway Genes).               | Ramón-Vidal, D. et al.*Database* Volume 2020, baaa090 (2020). [DOI: 10.1093/database/baaa090](https://doi.org/10.1093/database/baaa090)   |
| **iCW773 Metabolic Model** | Genome-scale metabolic network reconstruction of *C. glutamicum*.                                 | Zhang, Y. et al. *Biotechnology for Biofuels* 10, 169 (2017). [DOI: 10.1186/s13068-017-0856-3](https://doi.org/10.1186/s13068-017-0856-3)   |
| **iCGB21FR Model**         | Genome-scale metabolic model reconstruction (Support context).                                     | Feierabend, M. et al. *Frontiers in Microbiology* 12, 750206 (2021). [DOI: 10.3389/fmicb.2021.750206](https://doi.org/10.3389/fmicb.2021.750206) |
| **ecCGL1 Model**           | Enzyme-constrained metabolic model containing macromolecular parameters.                           | Niu, J. et al. *Biomolecules* 12, 1499 (2022). [DOI: 10.3390/biom12101499](https://doi.org/10.3390/biom12101499)                            |
| **PRODORIC**               | Position Weight Matrices (PWMs) database for prokaryotic promoter scanners.                        | Dudek, C.A. et al.*Nucleic Acids Research* 48, D1, D322–D329 (2020). [DOI: 10.1093/nar/gkz945](https://doi.org/10.1093/nar/gkz945)        |
| **STRING PPI**             | Protein-Protein physical and functional interaction network scores.                                | Szklarczyk, D. et al.*Nucleic Acids Research* 49, D1, D605–D612 (2021). [DOI: 10.1093/nar/gkaa1074](https://doi.org/10.1093/nar/gkaa1074) |
| **BRENDA**                 | Curated enzyme information system containing experimental kcat turnover numbers.                   | Chang, A. et al.*Nucleic Acids Research* 49, D1, D598–D604 (2021). [DOI: 10.1093/nar/gkaa1070](https://doi.org/10.1093/nar/gkaa1070)      |
| **iModulon DB**            | Independent component analysis (ICA) co-expression groups.                                         | Rychel, K. et al.*Nucleic Acids Research* 49, D1, D590–D597 (2021). [DOI: 10.1093/nar/gkaa1009](https://doi.org/10.1093/nar/gkaa1009)     |

---

### 4. Local Development & Setup

1. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
2. **Start the Local Server**:
   ```bash
   python run_server.py
   ```

   Open your browser and navigate to: `http://localhost:8000/index.html`.
3. **AI Summary API Configuration**:
   Configure API keys in the details panel of the right sidebar. Key parameters are stored locally in the browser's `localStorage`.

---

## 中文版本

## 谷氨酸棒状杆菌调控-代谢整合探索平台 (Cgl Regulation Explorer v0.5.0)

这是一个专门针对**谷氨酸棒状杆菌** (Corynebacterium glutamicum DSM 20300 / ATCC 13032) 研发的交互式、可复现的调控-代谢网路分析与仿真平台。

**线上正式访问地址**: [https://cgl-regulation.vercel.app/](https://cgl-regulation.vercel.app/)

该系统深度整合了转录因子 (TF)、sRNA 调控网络证据，与基因组级代谢反应模型、酶蛋白质约束参数，旨在为合成生物学工程改造提供靶点预测与设计方案。

---

### 1. 平台核心工作流

1. **调控证据整合 (Regulatory Evidence)**：汇集并标准化文献、启动子 Motifs、ChIP-seq 结合位点以及 sRNA-mRNA 互作数据。
2. **随机森林可信度评分 (RF Scoring)**：使用机器学习模型，对每一条调控关系进行打分。
3. **代谢网络映射 (Model Mapping)**：将靶基因关联至基因组级代谢模型 *iCW773* 中的酶、生化反应式及子通路。
4. **ecCGL1 酶约束定量 (Enzyme Annotations)**：补充代谢反应的动力学常数 kcat、分子量 (MW) 及 EC 编号，建立酶约束。
5. **通路级上游追溯 (Pathway-centered)**：自底向上，剖析控制特定代谢模块的全部潜在上游转录因子。
6. **工程候选靶点推荐 (Engineering Targets)**：计算“全局代谢影响力得分” (Global Metabolic Impact Score) 对 TF 进行排名。
7. **系统拓扑属性标注 (Topology Roles)**：结合 **Abasy Atlas** 识别全局枢纽基因 (Global Hubs)、模块化调控子 (Modular Regulators)。
8. **自动化案例文献报告 (AI Summary)**：通过 RAG 本地文献库与大语言模型，一键生成学术级的多角度功能综合报告。

---

### 2. 网站功能版块介绍

- **基因与转录因子浏览器 (Gene / TF Explorer)**：
  * 支持输入单基因或批量基因，构建级联调控子网络（Cytoscape.js 展示）。
  * 动态展示物理基因组定位、操纵子 (Operon) 结构及转录方向。
- **热力学启动子扫描 (Thermodynamic Promoter Scanner)**：
  * 支持交互式滑动窗口扫描，使用位置权重矩阵 (PWMs) 精确计算并图表化展示转录因子与 DNA 启动子区的热力学结合能分布。
- **代谢通路全局图谱 (Pathway View)**：
  * 可视化特定代谢通路的底物-产物反应流拓扑，追溯并按调控频数排列控制该通路的所有上游 TF。
- **在线代谢流与酶屏障仿真 (Active Simulation)**：
  * 支持在线运行动态调控通量平衡分析 (rFBA) 及酶约束通量平衡分析 (ecFBA)。
  * **酶瓶颈诊断 (Enzyme Bottlenecks)**：提供不同温度下（通过阿伦尼乌斯公式修正）的酶容量限制、阴影价格（Shadow Prices）分析，定位高风险代谢限速步骤。
- **网络拓扑分析 (Network Topology)**：
  * 自动计算当前子网络的度分布、介数中心性，并以三色温标动态警示系统表型级联风险（全局 Hub 为红色，模块化 TF 为蓝色，通路酶基因为绿色）。
- **iModulon 协同表达分析 (iModulon Explorer)**：
  * 引入独立成分分析 (ICA) 得到的 iModulon 基因群，展现基因在不同实验条件下的动态协同响应。
- **AI 智能综述与学术总结 (RAG Synthesis)**：
  * 集成了可定制 of AI 提供商设置（Gemini、OpenAI、DeepSeek 等），使用 RAG 搜索本地文献，提供关于靶基因的学术综述。

---

### 3. 数据来源及参考文献

在发表科研论文或使用本平台生成的图表、靶点方案时，请务必引用以下原始数据库及核心文献：

| 原始数据库 / 计算模型    | 数据类型与描述                                                   | 原始参考文献及链接                                                                                                                        |
| :----------------------- | :--------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------- |
| **CoryneRegNet 7** | 棒状杆菌属转录调控网络的核心骨架（转录因子、结合位点、靶基因）。 | Parise, M.T.D. et al.*Scientific Data* 7, 142 (2020). [DOI: 10.1038/s41597-020-0481-0](https://doi.org/10.1038/s41597-020-0481-0)        |
| **Abasy Atlas**    | 系统级调控拓扑角色映射（区分 Global Hub、Modular 等）。          | Ramón-Vidal, D. et al.*Database* (2020). [DOI: 10.1093/database/baaa090](https://doi.org/10.1093/database/baaa090)                      |
| **iCW773 模型**    | 谷氨酸棒状杆菌基因组级经典代谢反应模型。                         | Zhang, Y. et al. *Biotechnology for Biofuels* 10, 169 (2017). [DOI: 10.1186/s13068-017-0856-3](https://doi.org/10.1186/s13068-017-0856-3) |
| **iCGB21FR 模型**  | 谷氨酸棒状杆菌基因组级高解析重建代谢模型（辅助比对背景）。       | Feierabend, M. et al. *Frontiers in Microbiology* 12, 750206 (2021). [DOI: 10.3389/fmicb.2021.750206](https://doi.org/10.3389/fmicb.2021.750206) |
| **ecCGL1 模型**    | 结合了大分子拥挤和酶动力学约束（kcat/MW）的代谢模型。            | Niu, J. et al. *Biomolecules* 12, 1499 (2022). [DOI: 10.3390/biom12101499](https://doi.org/10.3390/biom12101499)                          |
| **PRODORIC**       | 原核生物转录因子结合位点及位置权重矩阵（PWMs）库。               | Dudek, C.A. et al.*Nucleic Acids Research* 48, D322–D329 (2020). [DOI: 10.1093/nar/gkz945](https://doi.org/10.1093/nar/gkz945)          |
| **STRING PPI**     | 蛋白质-蛋白质物理与功能相互作用可信度评分。                      | Szklarczyk, D. et al.*Nucleic Acids Research* 49, D605–D612 (2021). [DOI: 10.1093/nar/gkaa1074](https://doi.org/10.1093/nar/gkaa1074)   |
| **BRENDA**         | 酶学核心数据库，提供实验测定的代谢酶催化常数（kcat/Km）。        | Chang, A. et al.*Nucleic Acids Research* 49, D598–D604 (2021). [DOI: 10.1093/nar/gkaa1070](https://doi.org/10.1093/nar/gkaa1070)        |
| **iModulon DB**    | 基于ICA分解得到的转录协同表达模块数据库。                        | Rychel, K. et al.*Nucleic Acids Research* 49, D590–D597 (2021). [DOI: 10.1093/nar/gkaa1009](https://doi.org/10.1093/nar/gkaa1009)       |

---

### 4. 本地启动与部署步骤

1. **安装 Python 依赖**：
   ```bash
   pip install -r requirements.txt
   ```
2. **运行本地服务**：
   ```bash
   python run_server.py
   ```
   启动后，在浏览器访问：`http://localhost:8000/index.html` 即可进入控制台。
3. **AI 密钥配置**：
   直接在系统右侧详情面板的 “AI Provider Settings” 中填入个人的 API Key（仅保存在浏览器本地），点击 Refresh 即可生成报告。

---

### 5. Version History & Changelog / 版本历史与更新日志

#### v0.5.0 (July 2026 / 2026年7月)
* **English**:
  - Implemented formal GPR (Gene-Protein-Reaction) boolean rule evaluation (min/max operators) in dynamic rFBA simulations.
  - Exposed formerly silent thermodynamic directionality Safety Guard rollbacks as "Thermo-Stoichiometric Conflicts" in the UI.
  - Implemented mixed Chinese/English character and bigram tokenizer in RAG literature search.
  - Re-nested FBA model context in rFBA solver calls, speeding up simulation loops by 10-20x.
* **中文**:
  - 动态 rFBA 仿真全面引入 GPR (基因-蛋白质-反应) 精准布尔解析（min/max 算子），真实模拟多亚基复合酶的木桶限速效应和同工酶代偿。
  - 前端公开展示热力学单向约束冲突警告（被安全守护程序安全回滚的反应列表），辅助识别代谢模型结构缺陷。
  - RAG 文献检索引入双语混合分词器，完美支持中文及中英混合文献的高精相似度召回。
  - 重构 rFBA 通量计算引擎，实现单一上下文下的求解器复用，计算耗时缩短 10-20 倍。

#### v0.4.0 (July 2026 / 2026年7月)
* **English**:
  - Integrated **Abasy Atlas** systemic topology properties categorizing global hubs, modular regulators, and pathway genes.
  - Upgraded **RAG-based AI Literature summary** to natively parse Markdown Tables and render **Mermaid.js** flowcharts.
  - Added **Locus Link Anchoring** which turns raw tags in summaries into interactive links that center, zoom, and highlight Cytoscape nodes.
  - Linked search queries dynamically to highlight inputs and upstream regulators in the Global Metabolic Impact table.
  - Tightened dashboard layouts, sidebar spacing, and unified card typography/fonts.
* **中文**:
  - 深度整合 **Abasy Atlas** 拓扑角色分类，高亮展示全局枢纽、模块调节子及通路基因。
  - 升级 **RAG AI 智能综述**，原生支持 Markdown 表格转译与 **Mermaid.js** 交互式调控流程图。
  - 新增 **基因标签超级锚定**，AI 文本中的 locus tag 可一键点击联动 Cytoscape 画布缩放高亮及右侧详情页跳转。
  - 全局检索联动高亮，在 Global Metabolic Impact 列表中实时突出搜索靶基因的直接/上游转录因子。
  - 收紧侧边栏边距，移除冗余卡片并统一了全局学术字体排版。

#### v0.3.0 (June 2026 / 2026年6月)
* **English**: Upgraded simulator supporting DNA-binding energy curves (Chart.js), active ecFBA capacity bottleneck shadow price calculations, and interactive metabolic flow networks.
* **中文**: 升级启动子 DNA 结合能扫描曲线、ecFBA 酶瓶颈阻碍诊断与代谢通量分析流向图。
