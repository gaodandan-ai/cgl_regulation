# 实验室内网服务器部署指南 (Intranet Deployment Guide for 172.16.2.105)

本文档说明如何在课题组内网服务器 `172.16.2.105` 上部署与维护包含内部未发布数据（ChIP-seq 调控边 + 703 条件表达谱 + Peak 比对轨）的内网网站版本。

---

## 📋 部署环境要求 (Requirements)

- **服务器 IP**：`172.16.2.105`
- **监听端口**：`8010`（浏览器访问：`http://172.16.2.105:8010`）
- **推荐操作系统**：Ubuntu 20.04 / 22.04 LTS 或 Windows Server / Windows 10/11
- **推荐部署方式**：Docker & Docker Compose（或原生 Python 3.9+ 环境）

---

## 🚀 方式一：Linux / macOS 服务器一键部署 (推荐)

在 `172.16.2.105` 服务器终端中执行以下命令：

```bash
# 1. 一键部署并拉起服务
bash scripts/deploy/deploy_172.16.2.105.sh

# 2. 查看服务运行状态
bash scripts/deploy/deploy_172.16.2.105.sh status

# 3. 查看实时运行日志
bash scripts/deploy/deploy_172.16.2.105.sh logs

# 4. 停止服务
bash scripts/deploy/deploy_172.16.2.105.sh stop
```

---

## 🪟 方式二：Windows 服务器一键部署

在 `172.16.2.105` Windows 服务器上双击运行或在 CMD 中执行：

```cmd
scripts\deploy\deploy_172.16.2.105.bat
```

脚本会自动检测 Docker。如果已安装 Docker Desktop 则通过 Docker 容器运行；否则自动在控制台中通过 Python 启动 FastAPI 服务器。

---

## 🐳 方式三：Docker Compose 手动启动

如果您习惯直接使用 Docker 命令：

```bash
# 启动内网服务
docker compose -f scripts/deploy/docker-compose.intranet.yml up -d --build

# 检查健康状态
docker compose -f scripts/deploy/docker-compose.intranet.yml ps

# 查看日志
docker compose -f scripts/deploy/docker-compose.intranet.yml logs -f
```

---

## ⚙️ 内网配置说明 (Configuration)

在 `172.16.2.105` 上部署时，系统会自动应用以下环境变量：

| 环境变量 | 设定值 | 说明 |
| :--- | :--- | :--- |
| `PORT` | `8010` | 服务的默认监听端口 |
| `CGL_HOST` | `0.0.0.0` | 绑定所有网卡，允许局域网其他电脑访问 |
| `HEADLESS` | `true` | 服务器模式，不自动弹出浏览器 |
| `CGL_PUBLIC_DEPLOYMENT` | `false` | 保留全量 Cobra 代谢模型计算、本地 Ollama 大模型接入等功能 |
| `CGL_INTRANET_SERVER` | `172.16.2.105` | 激活网页 Header 紫金色 **`🧪 实验室内网版 (172.16.2.105)`** 徽章 |

---

## 🔒 数据更新与重新打包

当课题组在 `chip_rnaseq/` 目录下新增了 ChIP-seq 或 RNA-seq 样本数据时，只需在服务器上重新运行部署脚本：

```bash
bash deploy_172.16.2.105.sh
```

数据导入脚本会自动提取最新的多组学数据，重建 SQLite 数据库并热重载服务！
