# -*- mode: python ; coding: utf-8 -*-
import os

# ─── 基础数据文件 ─────────────────────────────────────────────────────────────
datas = [('web', 'web')]

if os.path.isdir('data/reference'):
    datas += Tree(
        'data/reference',
        prefix='data/reference',
        excludes=[
            '*.prebuild_*.bak',
            '*.db-wal',
            '*.db-shm',
            '*.db-journal',
        ],
    )

if os.path.isdir('backend/models'):
    datas.append(('backend/models', 'backend/models'))

# ─── 收集 pywebview 包数据（DLL、JS 资源等）─────────────────────────────────
from PyInstaller.utils.hooks import collect_all, collect_data_files

wv_datas, wv_bins, wv_hidden = collect_all('webview')
datas   += wv_datas
binaries = wv_bins   # webview 在 Windows 上通常包含 WebView2Loader.dll

# ─── Hidden Imports ────────────────────────────────────────────────────────────
HIDDEN = [
    # ── uvicorn ──
    'uvicorn.protocols.http.h11_impl',
    'uvicorn.loop.asyncio',
    'uvicorn.lifespan.on',
    'uvicorn.lifespan.off',
    'uvicorn.loop.auto',
    # ── scikit-learn ──
    'sklearn.ensemble._forest',
    'sklearn.utils._typedefs',
    # ── fastapi / pydantic ──
    'fastapi',
    'pydantic',
    # ── biology ──
    'cobra',
    'depinfo',
    # ── GUI ──
    'tkinter',
    'tkinter.ttk',
    # ── pywebview (Windows: pythonnet + WinForms WebView2 backend) ──
    'webview',
    'webview.platforms.winforms',   # Windows WebView2 via .NET WinForms
    'webview.platforms.edgechromium',
    'clr',                          # pythonnet CLR bridge
    'clr._extra',
    'System',
    'System.Windows.Forms',
    # ── backend modules ──
    'rag_service',
    'backend.app',
    'backend.gene_utils',
    'backend.kegg_client',
    'backend.metabolic_mapper',
    'backend.bio_handlers',
    'backend.sequence_tools',
    'backend.model_loader',
    'backend.thermo_pruner',
    'backend.simulation',
    'backend.schemas',
    'backend.objectives',
    'backend.thermodynamics',
    'backend.enzyme_thermal_params',
    # ── backend (flat-import aliases) ──
    'app',
    'gene_utils',
    'kegg_client',
    'metabolic_mapper',
    'bio_handlers',
    'sequence_tools',
    'model_loader',
    'thermo_pruner',
    'simulation',
    'schemas',
    'objectives',
    'thermodynamics',
    'enzyme_thermal_params',
]

HIDDEN += wv_hidden

# ─── Analysis ─────────────────────────────────────────────────────────────────
a = Analysis(
    ['launcher.pyw'],
    pathex=['backend'],
    binaries=binaries,
    datas=datas,
    hiddenimports=HIDDEN,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['equilibrator_api', 'tensorflow', 'torch'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='cgl_regulation',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,             # ← 无控制台窗口
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['icon.ico'],         # 任务栏/桌面图标
)
