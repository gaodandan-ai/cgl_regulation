语言: [ English ](README.md) | [ 中文版 ](README_zh.md)

# Cgl Regulation Explorer (v1.3.1)

**谷氨酸棒状杆菌系统生物学与合成生物学计算平台**

Cgl Regulation Explorer 是专为谷氨酸棒状杆菌（*Corynebacterium glutamicum* DSM 20300 / ATCC 13032）打造的综合计算平台，涵盖多维组学数据探索、调控网络分析、基因组规模代谢模型模拟以及 AI 辅助基因工程设计。

---

## 目录

- [系统概述](#系统概述)
- [核心计算与分析能力](#核心计算与分析能力)
  - [1. 多维组学调控网络可视化](#1-多维组学调控网络可视化)
  - [2. 5 轨高清互动基因组轨道浏览器](#2-5-轨高清互动基因组轨道浏览器)
  - [3. 环境依赖型组学与 iModulon 调控模块](#3-环境依赖型组学与-imodulon-调控模块)
  - [4. 基因组规模代谢模型与 FBA 通量模拟](#4-基因组规模代谢模型与-fba-通量模拟)
  - [5. 基于实测硬数据的 AI 代谢工程 Copilot](#5-基于实测硬数据的-ai-代谢工程-copilot)
- [代码库目录架构](#代码库目录架构)
- [安装与部署指南](#安装与部署指南)
  - [方式一：Windows 桌面客户端一键启动](#方式一windows-桌面客户端一键启动)
  - [方式二：开发服务器部署（Python FastAPI）](#方式二开发服务器部署python-fastapi)
  - [方式三：容器化部署（Docker）](#方式三容器化部署docker)
- [REST API 接口规范](#rest-api-接口规范)
- [主要组学数据来源](#主要组学数据来源)
- [测试与质量保障](#测试与质量保障)
- [引用与许可](#引用与许可)

---

## 系统概述

谷氨酸棒状杆菌（*Corynebacterium glutamicum*）是工业生物制造领域重要的氨基酸与有机酸生产菌株。Cgl Regulation Explorer 旨在建立原始转录调控数据集、基因组规模代谢网络与计算基因工程设计之间的桥梁。

平台整合了转录因子-靶基因（TF-TG）相互作用、sRNA-mRNA 转录后调控网络、操纵子结构、位置权重矩阵（PWM）、实测 ChIP-seq 结合峰、独立成分分析（ICA）iModulon 活泼度、STRING v12 蛋白质相互作用网络以及 iCW773 基因组规模代谢模型（GEM）。

---

## 核心计算与分析能力

### 1. 多维组学调控网络可视化
- **多层网络拓扑**：支持同时呈现转录因子、目标基因、sRNA 以及蛋白质相互作用（PPI）交叉网络。
- **视口纹理缓存**：启用 `textureOnViewport` 硬件加速，确保大规模网络图（>250 节点）拖拽顺畅。
- **Level-of-Detail (LOD) 分层渲染**：根据视口缩放比例动态调整文本标记可见性，在视口缩放与平移过程中保持 60 FPS 渲染帧率。

### 2. 5 轨高清互动基因组轨道浏览器
- **Track 1 (基因组坐标尺)**：呈现绝对碱基坐标与链方向刻度。
- **Track 2 (CDS 编码基因结构)**：以方向性 Block 符号呈现正链 (+) 与负链 (-) 编码基因结构。
- **Track 3 (Promoter 与 TSS 位点)**：标注转录起始位点（TSS）并高亮 Promoter 70bp 序列（包含 -35 与 -10 核心盒）。
- **Track 4 (ChIP-seq 结合峰信号强度)**：基于实测 ChIP-seq 与 CollectTF 结合得分绘制结合信号密度曲线。
- **Track 5 (sRNA/ncRNA 标记)**：展示非编码 RNA 基因位置。
- **交互控制**：支持视口放大、缩小、左右平移以及悬浮序列检视。

### 3. 环境依赖型组学与 iModulon 调控模块
- **iModulon 活泼度矩阵**：展示 9 种环境条件下 87 个 iModulon 活泼度矩阵及 F1-Score 覆盖对齐情况。
- **环境调控子网络**：支持筛选特定环境条件下的调控子网络（包括铁、氧、氮及逆境响应）。
- **图表双向交互**：支持图表区域框选并在主网络画布中实时高亮对应基因节点。

### 4. 基因组规模代谢模型与 FBA 通量模拟
- **iCW773 模型整合**：实现通量平衡分析（FBA）与代谢调整最小化（MOMA）算法。
- **基因扰动预测**：评估基因敲除与过表达对菌株生长速率与目标产物合成通量的影响。
- **热力学可行性检查**：结合自由能数据库对反应通量方向进行热力学约束校验。

### 5. 基于实测硬数据的 AI 代谢工程 Copilot
- **无幻觉组学上下文**：自动从 SQLite 数据库提取 RefSeq 坐标、Promoter 70bp 序列、TF 家族、效应物分子及 ChIP-seq 结合峰，用于硬数据约束提示词。
- **专业工程指令**：
  - `/design-crispri`：推荐用于 dCas9 基因敲降的 20bp gRNA 靶向区间与 PAM 位点（避开核心启动子区）。
  - `/solve-bottleneck`：诊断代谢瓶颈，推荐过表达与敲除协同靶点组合。
  - `/promoter-engineering`：基于效应物响应元件指导启动子突变文库构建。
- **多模型支持**：兼容 OpenAI、DeepSeek、Google Gemini 及本地离线 Ollama 模型。

---

## 代码库目录架构

```text
f:\cgl_regulation\
├── backend\                   # FastAPI 后端核心代码
│   ├── app.py                 # REST API 端点定义
│   ├── db_manager.py          # 线程安全的 SQLite 数据库管理器
│   ├── ai_handlers.py         # AI 组学硬事实与工程指令处理器
│   ├── bio_handlers.py        # 调控网络与通路分析
│   ├── graph_engine.py       # 网络拓扑与基序识别
│   ├── simulation.py          # FBA 与 MOMA 模拟算法
│   ├── model_loader.py        # SBML / iCW773 模型解析器
│   └── security.py            # 安全边界与请求校验
├── web\                       # 原生前端应用
│   ├── index.html             # 用户界面
│   ├── app.js                 # 前端应用逻辑
│   ├── style.css              # 样式定义
│   └── lib\                   # 模块化前端库
│       ├── genomicTrackBrowser.js  # 5 轨基因组轨道浏览器
│       ├── geneProfileViewer.js    # 360 度基因全景卡片
│       ├── networkTopology.js      # 网络拓扑模块
│       ├── icaConditionHeatmapView.js # iModulon 活泼度矩阵
│       ├── conditionRegulationView.js # 条件调控模块
│       ├── interventionPriorityView.js # 代谢工程靶点推荐
│       └── vendor\            # 本地第三方依赖库
├── data\                      # 数据库资产
│   ├── public_database.db     # SQLite 核心数据库
│   └── literature_cache.json  # RAG 文献缓存
├── data_pipeline\             # 数据 ETL 管道
│   ├── cli.py                 # 管道 CLI 入口
│   └── scripts\               # 组学解析脚本
├── tests\                     # Pytest API、数据、安全与模型测试
├── scripts\                   # 开发与维护脚本
│   ├── pipeline\              # GEO、ChIP-seq 与热力学数据构建
│   ├── analysis\              # 网络层级与跨株系分析
│   ├── build\                 # 构建与打包工具
│   └── archive\           # 归档历史脚本
├── launcher.pyw               # 桌面 GUI 启动器 (PyWebView)
├── run_server.py              # 开发服务器启动入口
├── cgl_regulation.spec        # PyInstaller 编译配置文件
├── requirements.txt           # Python 依赖清单
├── Dockerfile                 # 容器配置文件
└── README.md                  # 项目英文文档
```

---

## 安装与部署指南

### 环境要求
- **Python**: 3.10 或更高版本
- **操作系统**: Windows 10/11, macOS, Linux

---

### 方式一：跨平台桌面客户端一键启动 (Windows / macOS)

- **Windows 用户**：直接双击根目录下的启动脚本：
  ```cmd
  启动客户端.bat
  ```
- **macOS 用户**：在终端中运行 macOS 专属启动脚本：
  ```bash
  chmod +x start_mac.sh
  ./start_mac.sh
  ```
启动器将在后台开启 FastAPI 本地服务器并自动打开桌面容器（macOS 原生使用 WKWebView；Windows 使用 WebSockets/WebView2），全平台体验完全一致。

---

### 方式二：开发服务器部署（Python FastAPI）

1. 克隆代码库并安装依赖：
   ```bash
   git clone https://github.com/gaodandan-ai/cgl_regulation.git
   cd cgl_regulation
   pip install -r requirements.txt
   ```

2. 启动开发服务器：
   ```bash
   python run_server.py
   ```
   或使用 Uvicorn 启动：
   ```bash
   uvicorn backend.app:app --host 127.0.0.1 --port 8000 --reload
   ```

3. 访问界面：
   在浏览器中打开 `http://127.0.0.1:8000`。

---

### 方式三：容器化部署（Docker）

1. 构建 Docker 镜像：
   ```bash
   docker build -t cgl-regulation-explorer .
   ```

2. 运行容器：
   ```bash
   docker run -d -p 8000:8000 --name cgl-app cgl-regulation-explorer
   ```

---

## REST API 接口规范

| 端点 | 方法 | 说明 |
| :--- | :--- | :--- |
| `/api/gene/profile/{gene_id}` | `GET` | 查询目标基因 360 度全景 Profile 数据 |
| `/api/genomic_tracks/{gene_id}` | `GET` | 查询 5 轨基因组坐标与结合峰数据 |
| `/api/imodulon/condition` | `GET` | 查询指定环境下的 iModulon 活泼度矩阵 |
| `/api/intervention-targets` | `GET` | 返回代谢工程靶点推荐清单 |
| `/api/graph/motifs` | `GET` | 识别前馈环（FFL）网络基序 |
| `/api/ai/engineering_command` | `POST` | 执行 AI 代谢工程指令 (`/design-crispri`, `/solve-bottleneck`, `/promoter-engineering`) |
| `/api/check-update` | `GET` | 检查本地应用版本状态 |
| `/api/provenance` | `GET` | 查看数据库版本、数据哈希、发布清单及结果解释限制 |

---

## 主要组学数据来源

1. **CollectTF & RegPrecise**: 经过实验验证的转录因子结合位点与 PWM 矩阵。
2. **Abasy Atlas**: 系统级调控角色与多效性分类。
3. **iModulonDB**: 基于谷氨酸棒状杆菌 RNA-seq 数据集的独立成分分析集。
4. **STRING v12**: 蛋白质相互作用评分网络。
5. **NCBI RefSeq (NC_003450.3)**: 染色体基因组坐标与链方向标注。
6. **iCW773 GEM**: 谷氨酸棒状杆菌基因组规模代谢模型。

---

## 测试与质量保障

代码库包含覆盖后端 API 端点、数据库操作和计算逻辑的自动化测试套件：

运行测试套件：
```bash
python -m pytest
```

向 `main` 推送或提交 Pull Request 时，CI 会在 Python 3.10、3.11 和 3.13
运行完整测试。科研输出还可通过
`python scripts/analysis/validate_scientific_outputs.py` 执行机器可读的完整性审计。
数据证据分类、可复现元数据和验证边界详见
[数据来源与科研解释说明](docs/DATA_PROVENANCE.md)。

---

## 引用与许可

如果在科研出版物中使用了本计算平台，请进行如下引用：

```bibtex
@software{cgl_regulation_2026,
  author = {Gao, Dandan and DeepMind Agentic Coding Team},
  title = {Cgl Regulation Explorer: A Systems Biology and Synthetic Biology Platform for Corynebacterium glutamicum},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/gaodandan-ai/cgl_regulation}},
  version = {v1.3.1}
}
```

本软件遵循 **MIT License** 开源许可协议。详情请参阅 [LICENSE](LICENSE) 文件。
