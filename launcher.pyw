#!/usr/bin/env python3
"""
launcher.pyw
============
无控制台窗口的启动器：
  1. 显示美化的 tkinter Splash 加载界面
  2. 后台启动 FastAPI / Uvicorn 服务器
  3. 等待服务器就绪后，用 Chrome/Edge --app 模式打开（无地址栏）
  4. 服务器运行期间进程保持活跃；关闭浏览器后不关服务器（科研用途）

在 PyInstaller 打包时配置 console=False，icon=icon.ico
"""

import os, sys, time, socket, threading, subprocess, urllib.request

# ─── 路径 ────────────────────────────────────────────────────────────────────
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    ROOT_DIR = sys._MEIPASS
else:
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

PORT = int(os.environ.get("PORT", 8000))
URL  = f"http://127.0.0.1:{PORT}/index.html"

# ─── Splash 界面 ──────────────────────────────────────────────────────────────
def make_splash():
    """创建并返回一个无边框半透明的 Splash 窗口"""
    import tkinter as tk
    from tkinter import ttk

    root = tk.Tk()
    root.overrideredirect(True)          # 无标题栏、无边框
    root.configure(bg="#0f172a")
    root.attributes("-topmost", True)
    root.attributes("-alpha", 0.97)

    # 居中显示
    W, H = 480, 280
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    root.geometry(f"{W}x{H}+{(sw-W)//2}+{(sh-H)//2}")

    # ── Icon (icon.ico) ──────────────────────────────────────────────────────
    ico_path = os.path.join(ROOT_DIR, "icon.ico")
    if os.path.exists(ico_path):
        try:
            root.iconbitmap(ico_path)
        except Exception:
            pass

    # ── 标题标 ────────────────────────────────────────────────────────────────
    tk.Label(
        root,
        text="Cgl Regulation Explorer",
        font=("Segoe UI", 18, "bold"),
        fg="#e2e8f0",
        bg="#0f172a",
    ).pack(pady=(36, 0))

    tk.Label(
        root,
        text="Corynebacterium glutamicum Regulatory Network",
        font=("Segoe UI", 9),
        fg="#64748b",
        bg="#0f172a",
    ).pack(pady=(4, 0))

    # ── 版本 ─────────────────────────────────────────────────────────────────
    try:
        import json
        vf = os.path.join(ROOT_DIR, "web", "version.json")
        with open(vf, encoding="utf-8") as f:
            ver = json.load(f).get("version", "")
        ver_text = f"v{ver}"
    except Exception:
        ver_text = ""

    if ver_text:
        tk.Label(
            root,
            text=ver_text,
            font=("Segoe UI", 8),
            fg="#334155",
            bg="#0f172a",
        ).pack()

    # ── 进度条 ────────────────────────────────────────────────────────────────
    frame = tk.Frame(root, bg="#0f172a")
    frame.pack(fill="x", padx=48, pady=(32, 0))

    style = ttk.Style()
    style.theme_use("clam")
    style.configure(
        "Cgl.Horizontal.TProgressbar",
        troughcolor="#1e293b",
        background="#0ea5e9",
        darkcolor="#0ea5e9",
        lightcolor="#38bdf8",
        bordercolor="#0f172a",
        thickness=6,
    )

    pbar = ttk.Progressbar(
        frame,
        style="Cgl.Horizontal.TProgressbar",
        orient="horizontal",
        length=384,
        mode="indeterminate",
    )
    pbar.pack()
    pbar.start(12)

    # ── 状态文字 ──────────────────────────────────────────────────────────────
    status_var = tk.StringVar(value="正在启动服务器...")
    tk.Label(
        root,
        textvariable=status_var,
        font=("Segoe UI", 9),
        fg="#475569",
        bg="#0f172a",
    ).pack(pady=(12, 0))

    # ── 版权 ─────────────────────────────────────────────────────────────────
    tk.Label(
        root,
        text="© 2026 gaodandan-ai",
        font=("Segoe UI", 7),
        fg="#1e293b",
        bg="#0f172a",
    ).pack(side="bottom", pady=12)

    return root, pbar, status_var


# ─── 服务器可达性检测 ─────────────────────────────────────────────────────────
def wait_for_server(host="127.0.0.1", port=PORT, timeout=60):
    """轮询直到服务器端口可用"""
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                # 额外验证 HTTP 端点响应
                try:
                    urllib.request.urlopen(f"http://{host}:{port}/", timeout=1)
                    return True
                except Exception:
                    return True   # 端口通了就算 OK
        except OSError:
            time.sleep(0.25)
    return False


# ─── 浏览器（App 模式，无地址栏）─────────────────────────────────────────────
CHROME_PATHS = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
]
EDGE_PATHS = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]

def find_browser():
    for p in CHROME_PATHS + EDGE_PATHS:
        if os.path.exists(p):
            return p
    return None

def open_app_window(url, browser_path, icon_path):
    """
    --app=URL   : 无地址栏、无标签页条的 App 窗口
    --window-size / --window-position: 设置初始大小
    --user-data-dir: 独立 profile，避免和普通 Chrome 冲突
    """
    profile_dir = os.path.join(os.path.expanduser("~"), ".cgl_regulation_browser_profile")
    args = [
        browser_path,
        f"--app={url}",
        f"--user-data-dir={profile_dir}",
        "--window-size=1400,900",
        "--window-position=60,40",
        "--disable-extensions",
        "--no-first-run",
        "--disable-default-apps",
    ]
    subprocess.Popen(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# ─── 后台启动 uvicorn ─────────────────────────────────────────────────────────
def start_server_background():
    """在后台线程里启动 uvicorn；不阻塞 GUI 主线程"""
    import uvicorn
    # Redirect stdout/stderr to avoid any console popup
    sys.stdout = open(os.devnull, "w")
    sys.stderr = open(os.devnull, "w")
    try:
        from backend.app import app as fastapi_app
        uvicorn.run(fastapi_app, host="0.0.0.0", port=PORT, log_level="warning")
    except Exception:
        pass


# ─── 主流程 ───────────────────────────────────────────────────────────────────
def main():
    import tkinter as tk

    # 1. 先检查端口是否已经在用（避免重复启动）
    already_running = False
    try:
        with socket.create_connection(("127.0.0.1", PORT), timeout=0.3):
            already_running = True
    except OSError:
        pass

    # 2. 显示 Splash
    root, pbar, status_var = make_splash()

    def run():
        nonlocal already_running

        if not already_running:
            # 3. 后台启动服务器
            server_thread = threading.Thread(target=start_server_background, daemon=True)
            server_thread.start()

            status_var.set("正在加载数据模型...")
            ok = wait_for_server(timeout=90)
            if not ok:
                status_var.set("启动超时，请检查环境配置")
                time.sleep(3)
                root.quit()
                return
        else:
            status_var.set("检测到已在运行，直接打开...")
            time.sleep(0.5)

        status_var.set("正在打开界面...")

        # 4. App 模式浏览器
        browser = find_browser()
        icon_path = os.path.join(ROOT_DIR, "icon.ico")
        if browser:
            open_app_window(URL, browser, icon_path)
        else:
            import webbrowser
            webbrowser.open(URL)

        time.sleep(1.2)
        # 5. 关闭 Splash
        root.quit()

    threading.Thread(target=run, daemon=True).start()
    root.mainloop()

    try:
        root.destroy()
    except Exception:
        pass

    # 6. 保持主进程活跃（服务器在后台线程里跑）
    # 服务器线程是 daemon，所以这里 join 主线程无限等待
    try:
        while True:
            time.sleep(60)
    except (KeyboardInterrupt, SystemExit):
        pass


if __name__ == "__main__":
    main()
