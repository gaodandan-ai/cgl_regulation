#!/usr/bin/env python3
"""
launcher.pyw
============
无控制台窗口的启动器 (v4 — pywebview 原生窗口):
  1. Windows 命名互斥锁确保只运行一个实例
  2. 检测端口是否已被占用 → 直接复用，不重新绑定
  3. 后台启动 FastAPI / Uvicorn 服务器（仅当端口空闲时）
  4. 显示美化的 tkinter Splash 加载界面（服务器就绪前）
  5. 服务器就绪后切换到 pywebview 原生窗口（不依赖任何外部浏览器）
  6. PID 文件 + atexit + 信号处理，确保进程彻底退出
  7. pywebview 窗口关闭后自动关闭服务器
"""

import os, sys, time, socket, threading, urllib.request, atexit, signal, json, logging, traceback

# ─── 路径 ────────────────────────────────────────────────────────────────────
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    ROOT_DIR = sys._MEIPASS
else:
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

PORT = int(os.environ.get("PORT", 8000))
_app_data_root = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
LOG_DIR = os.path.join(_app_data_root, "Cgl Regulation Explorer", "logs")
LOG_FILE = os.path.join(LOG_DIR, "launcher.log")

# ─── PID 文件路径（用于优雅退出 & 安装程序检测）────────────────────────────────
PID_FILE = os.path.join(os.path.expanduser("~"), ".cgl_regulation_server.pid")


# ─── 0. 优雅退出基础设施 ─────────────────────────────────────────────────────

def _write_pid_file():
    """将当前 PID 写入文件，供其他进程或安装程序识别"""
    try:
        with open(PID_FILE, "w") as f:
            f.write(str(os.getpid()))
    except Exception:
        pass


def _remove_pid_file():
    """退出时清理 PID 文件"""
    try:
        if os.path.exists(PID_FILE):
            os.remove(PID_FILE)
    except Exception:
        pass


def _graceful_shutdown(signum=None, frame=None):
    """统一退出入口：清理资源后强制退出，防止 daemon 线程阻塞"""
    _remove_pid_file()
    os._exit(0)


def _setup_exit_hooks():
    """注册 atexit 和信号处理器"""
    atexit.register(_remove_pid_file)
    try:
        signal.signal(signal.SIGTERM, _graceful_shutdown)   # kill 命令
    except (OSError, ValueError):
        pass
    try:
        signal.signal(signal.SIGBREAK, _graceful_shutdown)  # Windows Ctrl+Break
    except (OSError, AttributeError):
        pass


# ─── 1. 单实例检测（Windows Named Mutex）─────────────────────────────────────
def acquire_single_instance_lock():
    """
    尝试创建一个 Windows 命名互斥锁。
    如果锁已存在，说明另一个实例正在运行 → 只打开新窗口，然后退出。
    返回 (mutex_handle_or_None, is_first_instance)
    """
    try:
        import ctypes
        kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        mutex = kernel32.CreateMutexW(None, True, "CglRegulationExplorer_SingleInstance_Mutex")
        last_err = ctypes.get_last_error()
        ERROR_ALREADY_EXISTS = 183
        if last_err == ERROR_ALREADY_EXISTS:
            return mutex, False   # 已有实例在运行
        return mutex, True        # 第一个实例
    except Exception:
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
    """Verify the application identity instead of trusting port occupancy."""
    try:
        req = urllib.request.Request(
            f"http://127.0.0.1:{port}/api/health",
            headers={"User-Agent": "CglLauncher/4"},
        )
        with urllib.request.urlopen(req, timeout=1) as r:
            payload = json.loads(r.read().decode("utf-8"))
            return r.status == 200 and payload.get("app") == "cgl-regulation"
    except Exception:
        return False


# ─── 3. Splash 界面（tkinter，仅在服务器就绪前显示）────────────────────────────
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
    """Start the loopback-only API and retain diagnostics in a rotating log."""
    try:
        import uvicorn
        from logging.handlers import RotatingFileHandler

        os.makedirs(LOG_DIR, exist_ok=True)
        handler = RotatingFileHandler(LOG_FILE, maxBytes=2_000_000, backupCount=3, encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s"))
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.INFO)
        root_logger.addHandler(handler)

        backend_dir = os.path.join(ROOT_DIR, "backend")
        if backend_dir not in sys.path:
            sys.path.insert(0, backend_dir)
        if ROOT_DIR not in sys.path:
            sys.path.insert(0, ROOT_DIR)

        from backend.app import app as fastapi_app
        uvicorn.run(fastapi_app, host="127.0.0.1", port=port, log_level="info")
    except Exception:
        try:
            os.makedirs(LOG_DIR, exist_ok=True)
            with open(LOG_FILE, "a", encoding="utf-8") as stream:
                stream.write(traceback.format_exc())
        except Exception:
            pass


def wait_for_server(port: int, timeout: int = 120) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        if is_our_server_ready(port):
            return True
        time.sleep(0.3)
    return False


# ─── 5. pywebview 原生窗口 ────────────────────────────────────────────────────
def _app_title() -> str:
    """读取 version.json 构造窗口标题"""
    try:
        import json
        vf = os.path.join(ROOT_DIR, "web", "version.json")
        with open(vf, encoding="utf-8") as f:
            ver = json.load(f).get("version", "")
        if ver:
            return f"Cgl Regulation Explorer  v{ver}"
    except Exception:
        pass
    return "Cgl Regulation Explorer"


def open_native_window(url: str):
    """
    在主线程打开 pywebview 原生窗口，阻塞直到用户关闭。
    关闭后调用 _graceful_shutdown() 退出整个进程。
    """
    import webview

    ico_path = os.path.join(ROOT_DIR, "icon.ico")
    storage_path = os.path.join(os.path.expanduser("~"), ".cgl_regulation_webview")

    webview.create_window(
        _app_title(),
        url,
        width=1400,
        height=900,
        x=60,
        y=40,
        min_size=(900, 600),
        background_color="#0f172a",  # 与 Splash 背景一致，消除白闪
        confirm_close=False,
        text_select=True,            # 允许文本选择（论文、基因名等）
    )

    # private_mode=False → 保持 localStorage / cookie 跨会话
    # storage_path        → WebView2 数据存放位置（cookie、缓存等）
    # icon                → 任务栏图标（Windows .ico）
    webview.start(
        debug=False,
        private_mode=False,
        storage_path=storage_path,
        icon=ico_path if os.path.exists(ico_path) else None,
    )
    # webview.start() 在窗口关闭时返回 → 触发服务器退出
    _graceful_shutdown()


# ─── 主流程 ───────────────────────────────────────────────────────────────────
def main():
    # ── 优雅退出基础设施初始化 ────────────────────────────────────────────────
    _setup_exit_hooks()
    _write_pid_file()

    # ── 单实例检测 ────────────────────────────────────────────────────────────
    mutex, is_first = acquire_single_instance_lock()

    url = f"http://127.0.0.1:{PORT}/index.html"

    if not is_first:
        # 已有实例在运行 → 复用现有服务器，打开一个新 pywebview 窗口
        _remove_pid_file()
        try:
            import webview
            storage_path = os.path.join(os.path.expanduser("~"), ".cgl_regulation_webview")
            webview.create_window(
                _app_title(), url,
                width=1400, height=900,
                min_size=(900, 600),
                background_color="#0f172a",
            )
            webview.start(debug=False, private_mode=False, storage_path=storage_path)
        except Exception:
            import webbrowser
            webbrowser.open(url)
        return

    # ── 端口占用检测 ──────────────────────────────────────────────────────────
    port_busy = is_port_in_use(PORT)

    # ── Phase 1: tkinter Splash（主线程 mainloop）────────────────────────────
    root, pbar, status_var = make_splash()

    def _run_in_bg():
        """后台线程：启动服务器并等待就绪，就绪后退出 tkinter mainloop"""
        if port_busy:
            status_var.set(f"检测到服务已在运行（端口 {PORT}），正在连接...")
            ok = wait_for_server(PORT, timeout=10)
            if not ok:
                status_var.set(f"端口 {PORT} 被其他程序占用，请先关闭后重试")
                time.sleep(4)
                root.quit()
                return
        else:
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

        status_var.set("正在初始化界面...")
        time.sleep(0.4)
        root.quit()   # ← 结束 tkinter mainloop，主线程继续向下执行

    threading.Thread(target=_run_in_bg, daemon=True).start()
    root.mainloop()   # 主线程阻塞于此直到 root.quit() 被调用

    try:
        root.destroy()
    except Exception:
        pass

    # ── Phase 2: pywebview 原生窗口（主线程接管）────────────────────────────
    if is_our_server_ready(PORT):
        open_native_window(url)   # 阻塞到窗口关闭 → 内部调用 _graceful_shutdown()

    # ── Mutex 释放 ────────────────────────────────────────────────────────────
    if mutex:
        try:
            import ctypes
            ctypes.WinDLL('kernel32').CloseHandle(mutex)
        except Exception:
            pass

    _graceful_shutdown()


if __name__ == "__main__":
    main()
