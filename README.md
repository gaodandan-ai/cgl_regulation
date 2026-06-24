# Corynebacterium glutamicum Regulatory Network Explorer (C.g. RNE)

[简体中文](#简体中文) | [English](#english)

---

## 简体中文

谷氨酸棒状杆菌 (*Corynebacterium glutamicum*) 基因调控与 sRNA 调控网络交互式可视化分析与预测平台。

本系统集成了本地基因调控数据（TF-TG 关系、操纵子结构、预测的 sRNA 相互作用关系）与基于 Gemini 模型的 AI 智能学术分析助手，提供从拓扑可视化、表达扰动效应模拟到 AI 文献综述的一站式分析体验。

### 主要功能

1. **交互式网页端网络浏览器 (Web Explorer)**
   - **多模式输入**：支持逐个基因添加，或通过空格、换行、分号批量导入多个基因/sRNA 开展联合分析。
   - **精细调控过滤**：可动态过滤激活 (+)、抑制 (-)、双重/Sigma 调控以及 sRNA 预测调控，支持调节 CopraRNA Rank 限制阈值。
   - **关键拓扑分析**：一键筛选“只显示共同调控靶标（Co-regulated）”或“只显示转录因子（TF）靶标”，自动识别当前网络的核心 Hub 节点。
   - **高级画布交互**：基于 Cytoscape.js，提供 CoSE、同心圆、层级等多种布局，支持画布内基因搜索定位、高亮一阶调控子网，以及高分辨率 PNG / 调控明细 CSV 导出。

2. **基因表达扰动预测模拟 (Perturbation Simulator)**
   - 模拟对当前转录因子/sRNA 进行 **基因过表达 (OE)** 或 **基因敲除 (KO)**。
   - 动态预测并渲染下游靶基因的表达水平响应（表达增强 ⬆、表达减弱 ⬇、复杂/双重 ↕），并支持一键导出扰动预测 CSV 报表。

3. **AI 智能学术分析助手 (AI Assistant)**
   - **文献与功能一键总结**：自动检索 PubMed 数据库获取当前基因的最新文献摘要，并结合 Gemini 模型自动生成排版精美的学术综述。
   - **AI 基因分析助手**：输入功能或生理描述（如“抗逆相关调控因子”），AI 自动匹配并绘制关联调控网络。
   - **AI 通路分析助手**：输入代谢通路（如“赖氨酸合成”、“糖酵解”），AI 自动提取关键基因并一键生成调控子网。

4. **命令行可视化工具 (CLI Visualizer)**
   - 提供 Python 脚本 `visualize_network.py`。
   - 可一键生成符合学术期刊发表标准的高分辨率静态网络图（基于 Matplotlib / NetworkX）。
   - 可选生成基于 PyVis 的三维力导向交互式 HTML 网络图（支持节点悬浮属性详情显示）。

---

### 项目结构

```text
├── data/                          # 本地调控与基因数据库 (CSV/XLSX 格式)
│   ├── regulations.csv            # TF-TG 调控数据
│   ├── rna_regulation.csv         # sRNA-mRNA 相互作用数据
│   ├── operons.csv                # 操纵子结构数据
│   ├── gene_mapping.csv           # 基因名与 Locus Tag 映射关系表
│   └── gene_mapping.xlsx          # 原始 Excel 基因映射表
├── web/                           # 前端网页 Explorer 资源目录
│   ├── index.html                 # 网页端主页面
│   ├── style.css                  # 自定义网页样式
│   ├── app.js                     # 网页端交互与 API 请求核心逻辑
│   └── lib/                       # 前端第三方依赖库 (Cytoscape.js, PapaParse, TomSelect 等)
├── run_server.py                  # Python 本地轻量级 API 代理及文件服务器
├── visualize_network.py           # 命令行绘图与网络生成工具 (支持输出 PNG 及交互 HTML)
├── requirements.txt               # Python 依赖包声明
├── run.bat                        # Windows 一键配置启动脚本
└── run.sh                         # macOS/Linux 一键配置启动脚本
```

---

### 快速开始

#### 方法 A：使用一键启动脚本（推荐，开箱即用）

1. 双击运行项目根目录下的启动脚本：
   - **Windows**：双击 `run.bat`
   - **macOS / Linux**：在终端运行 `./run.sh`（首次运行可能需要执行 `chmod +x run.sh` 赋予权限）
2. 脚本会自动检测 Python 环境，安装所需依赖，并启动本地服务器，最后自动在默认浏览器中打开系统主页：
   [http://localhost:8000/index.html](http://localhost:8000/index.html)

#### 方法 B：手动运行

如果您希望手动配置和运行系统，请按照以下步骤操作：

1. **安装 Python 3 依赖**
   打开终端或 PowerShell，运行：
   ```bash
   pip install -r requirements.txt
   ```
2. **启动本地服务器**
   运行后台服务端脚本（处理本地静态页面托管及 PubMed/Gemini API 代理转发）：
   ```bash
   python run_server.py
   ```
3. **打开浏览器**
   在浏览器中访问：[http://localhost:8000/index.html](http://localhost:8000/index.html)

---

### 使用 AI 分析功能

为了启用 AI 学术总结、通路分析和基因提取助手，请按以下步骤配置您的 API Key：
1. 访问 [Google AI Studio](https://aistudio.google.com/app/apikey) 获取免费的 **Gemini API Key**。
2. 打开本系统网页，在右侧详情面板的 **AI 文献与功能智能总结** 区域粘贴您的 API Key，并点击 **保存**。
3. 密钥将安全地保存在您本地浏览器的 `localStorage` 中，不会上传至任何第三方服务器。

---

### 使用命令行可视化工具

除了网页端，您还可以直接使用 `visualize_network.py` 从终端生成高品质的调控网络图：

```bash
# 示例：绘制以 whiB4 为中心的一阶调控网络，同时输出 PNG 和 交互式 HTML
python visualize_network.py --gene whiB4 --steps 1

# 示例：仅展示受到 2 个或以上 TF 共同调控的靶基因
python visualize_network.py --gene sigH --steps 1 --only-coregulated

# 查看所有可用参数
python visualize_network.py --help
```

---
---

## English

An interactive visualization, simulation, and analysis platform for the transcriptional and sRNA-mediated regulatory network of *Corynebacterium glutamicum* DSM 20300 = ATCC 13032.

This system integrates local database relations (TF-TG links, operon structures, predicted sRNA interactions) with Gemini-powered AI assistants to provide a seamless research workflow from network topology exploration and genetic perturbation simulation to literature synthesis.

### Key Features

1. **Interactive Web-based Explorer**
   - **Multi-mode Querying**: Input a single gene or batch paste multiple genes/sRNAs (delimited by spaces, newlines, commas, or semicolons) to query.
   - **Granular Filters**: Dynamically toggle Activation (+), Repression (-), Dual/Sigma, and sRNA links. Fine-tune sRNA inputs using the CopraRNA rank threshold slider.
   - **Network Diagnostics**: Highlight "Only Co-regulated Targets" or "Only Transcription Factor Targets" instantly. Auto-calculate network hubs.
   - **Rich UI Controls**: Supporting layouts like CoSE, concentric, circle, and hierarchy. Includes in-canvas node search/flashing, 1st-degree neighbor highlighting, high-res PNG export, and CSV download.

2. **Perturbation Simulator**
   - Simulate **Overexpression (OE)** or **Knockout (KO)** of a TF or sRNA.
   - Graphically predict down-stream target responses (Upregulated ⬆, Downregulated ⬇, or Dual/Complex ↕) and export results as a CSV sheet.

3. **AI Research Assistant**
   - **Lit Synthesis**: Queries PubMed API for abstract data of the selected gene, combining it with Gemini API to compile an academic summary.
   - **AI Gene/Pathway Helper**: Input natural language (e.g., "Lysine synthesis pathway" or "stress response regulators"), and the AI will auto-extract locus tags and construct the corresponding subnetwork in one click.

4. **CLI Visualizer**
   - A command-line tool `visualize_network.py` built on Matplotlib and NetworkX.
   - Generates publication-ready static PNGs and optional interactive PyVis HTML files with floating tooltips.

---

### Quick Start

#### Method A: Click-and-Run Scripts (Recommended)

1. Double-click the startup script in the root directory:
   - **Windows**: `run.bat`
   - **macOS / Linux**: Run `./run.sh` in your terminal (make it executable first via `chmod +x run.sh` if needed).
2. The script will automatically check Python 3, install necessary libraries from `requirements.txt`, launch the server, and open the tool page in your default web browser:
   [http://localhost:8000/index.html](http://localhost:8000/index.html)

#### Method B: Manual Startup

1. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   ```
2. **Launch Server**
   ```bash
   python run_server.py
   ```
3. **Open Explorer**
   Open your browser and navigate to: [http://localhost:8000/index.html](http://localhost:8000/index.html)

---

### Activating AI Assistant Features

1. Obtain a free **Gemini API Key** from [Google AI Studio](https://aistudio.google.com/app/apikey).
2. In the Web Explorer UI, open the right sidebar panel, paste the key in the **AI Literature & Function Summary** section, and click **Save**.
3. Keys are stored safely in your browser's local storage and are sent directly to the local server proxy to request Gemini APIs.
