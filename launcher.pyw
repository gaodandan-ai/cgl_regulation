#!/usr/bin/env python3
"""
launcher.pyw
============
无控制台窗口的启动器 (v2 — 单实例 + 端口复用):
  1. Windows 命名互斥锁确保只运行一个实例
  2. 检测端口是否已被占用 → 直接复用，不重新绑定
  3. 后台启动 FastAPI / Uvicorn 服务器（仅当端口空闲时）
  4. 显示美化的 tkinter Splash 加载界面
  5. 服务器就绪后用 Chrome/Edge --app 模式打开（无地址栏）
"""

import os, sys, time, socket, threading, subprocess, urllib.request

# ─── 路径 ────────────────────────────────────────────────────────────────────
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    ROOT_DIR = sys._MEIPASS
else:
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

PORT = int(os.environ.get("PORT", 8000))


# ─── 1. 单实例检测（Windows Named Mutex）─────────────────────────────────────
def acquire_single_instance_lock():
    """
    尝试创建一个 Windows 命名互斥锁。
    如果锁已存在，说明另一个实例正在运行 → 只打开浏览器，然后退出。
    返回 (mutex_handle_or_None, is_first_instance)
    """
    try:
        import ctypes
        from ctypes import wintypes
        kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        mutex = kernel32.CreateMutexW(None, True, "CglRegulationExplorer_SingleInstance_Mutex")
        last_err = ctypes.get_last_error()
        ERROR_ALREADY_EXISTS = 183
        if last_err == ERROR_ALREADY_EXISTS:
            return mutex, False   # 已有实例在运行
        return mutex, True        # 第一个实例
    except Exception:
        # 无法创建 Mutex（非 Windows 环境）→ 允许启动
        return None, True


# ─── 2. 端口占用检测 ──────────────────────────────────────────────────────────
def is_port_in_use(port: int) -> bool:
    """检查 127.0.0.1:port 是否已有进程监听"""
    try:
        with socket.create_connection(("127.0.0.1", port), timeout=0.3):
            return True
    except OSError:
        return False


def is_our_server_ready(port: int) -> bool:
    """检查端口上运行的是否是我们的 FastAPI 服务（能返回 HTTP 响应）"""
    try:
        req = urllib.request.Request(
            f"http://127.0.0.1:{port}/",
            headers={"User-Agent": "CglLauncher/2"},
        )
        with urllib.request.urlopen(req, timeout=1) as r:
            return r.status < 500
    except Exception:
        return is_port_in_use(port)


# ─── 3. Splash 界面 ───────────────────────────────────────────────────────────
def make_splash():
    import tkinter as tk
    from tkinter import ttk

    root = tk.Tk()
    root.overrideredirect(True)
    root.configure(bg="#0f172a")
    root.attributes("-topmost", True)
    root.attributes("-alpha", 0.97)

    W, H = 480, 280
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    root.geometry(f"{W}x{H}+{(sw-W)//2}+{(sh-H)//2}")

    ico_path = os.path.join(ROOT_DIR, "icon.ico")
    if os.path.exists(ico_path):
        try:
            root.iconbitmap(ico_path)
        except Exception:
            pass

    tk.Label(
        root, text="Cgl Regulation Explorer",
        font=("Segoe UI", 18, "bold"), fg="#e2e8f0", bg="#0f172a",
    ).pack(pady=(36, 0))

    tk.Label(
        root, text="Corynebacterium glutamicum Regulatory Network",
        font=("Segoe UI", 9), fg="#64748b", bg="#0f172a",
    ).pack(pady=(4, 0))

    # 版本号
    try:
        import json
        vf = os.path.join(ROOT_DIR, "web", "version.json")
        with open(vf, encoding="utf-8") as f:
            ver = json.load(f).get("version", "")
        if ver:
            tk.Label(root, text=f"v{ver}", font=("Segoe UI", 8),
                     fg="#334155", bg="#0f172a").pack()
    except Exception:
        pass

    # 进度条
    frame = tk.Frame(root, bg="#0f172a")
    frame.pack(fill="x", padx=48, pady=(32, 0))

    style = ttk.Style()
    style.theme_use("clam")
    style.configure(
        "Cgl.Horizontal.TProgressbar",
        troughcolor="#1e293b", background="#0ea5e9",
        darkcolor="#0ea5e9", lightcolor="#38bdf8",
        bordercolor="#0f172a", thickness=6,
    )
    pbar = ttk.Progressbar(
        frame, style="Cgl.Horizontal.TProgressbar",
        orient="horizontal", length=384, mode="indeterminate",
    )
    pbar.pack()
    pbar.start(12)

    status_var = tk.StringVar(value="正在启动...")
    tk.Label(root, textvariable=status_var, font=("Segoe UI", 9),
             fg="#475569", bg="#0f172a").pack(pady=(12, 0))

    tk.Label(root, text="© 2026 gaodandan-ai", font=("Segoe UI", 7),
             fg="#1e293b", bg="#0f172a").pack(side="bottom", pady=12)

    return root, pbar, status_var


# ─── 4. 后台服务器 ─────────────────────────────────────────────────────────────
def start_server_background(port: int):
    """在后台线程启动 uvicorn；丢弃所有输出"""
    try:
        import uvicorn
        # 重定向输出避免弹出任何窗口
        null = open(os.devnull, "w")
        sys.stdout = null
        sys.stderr = null

        # 确保 backend 可 import
        backend_dir = os.path.join(ROOT_DIR, "backend")
        if backend_dir not in sys.path:
            sys.path.insert(0, backend_dir)
        if ROOT_DIR not in sys.path:
            sys.path.insert(0, ROOT_DIR)

        from backend.app import app as fastapi_app
        uvicorn.run(fastapi_app, host="0.0.0.0", port=port, log_level="error")
    except Exception:
        pass


def wait_for_server(port: int, timeout: int = 90) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        if is_our_server_ready(port):
            return True
        time.sleep(0.3)
    return False


# ─── 5. 浏览器（App 模式）────────────────────────────────────────────────────
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

def open_app_window(url: str, browser_path: str):
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


# ─── 主流程 ───────────────────────────────────────────────────────────────────
def main():
    # ── 单实例检测 ────────────────────────────────────────────────────────────
    mutex, is_first = acquire_single_instance_lock()

    url = f"http://127.0.0.1:{PORT}/index.html"

    if not is_first:
        # 另一个实例已在运行 → 只打开浏览器，直接退出
        browser = find_browser()
        if browser:
            open_app_window(url, browser)
        else:
            import webbrowser
            webbrowser.open(url)
        return

    # ── 端口占用检测 ──────────────────────────────────────────────────────────
    port_busy = is_port_in_use(PORT)

    # ── Splash ────────────────────────────────────────────────────────────────
    root, pbar, status_var = make_splash()

    def run():
        if port_busy:
            # 端口已被占用 → 复用现有服务器
            status_var.set(f"检测到服务已在运行（端口 {PORT}），正在连接...")
            # 稍等确保对方完全就绪
            ok = wait_for_server(PORT, timeout=10)
            if not ok:
                status_var.set(f"端口 {PORT} 被其他程序占用，请先关闭后重试")
                time.sleep(4)
                root.quit()
                return
        else:
            # 端口空闲 → 正常启动服务器
            status_var.set("正在启动服务器...")
            t = threading.Thread(
                target=start_server_background,
                args=(PORT,),
                daemon=True,
            )
            t.start()
            status_var.set("正在加载调控网络数据...")
            ok = wait_for_server(PORT, timeout=120)
            if not ok:
                status_var.set("启动超时，请检查 Python 环境")
                time.sleep(4)
                root.quit()
                return

        status_var.set("正在打开界面...")
        browser = find_browser()
        if browser:
            open_app_window(url, browser)
        else:
            import webbrowser
            webbrowser.open(url)

        time.sleep(1.0)
        root.quit()

    threading.Thread(target=run, daemon=True).start()
    root.mainloop()

    try:
        root.destroy()
    except Exception:
        pass

    # 保持进程活跃（后台服务器 daemon 线程需要主线程存活）
    if not port_busy:
        try:
            while True:
                time.sleep(60)
        except (KeyboardInterrupt, SystemExit):
            pass

    # 释放 Mutex
    if mutex:
        try:
            import ctypes
            ctypes.WinDLL('kernel32').CloseHandle(mutex)
        except Exception:
            pass


if __name__ == "__main__":
    main()
