# -*- mode: python ; coding: utf-8 -*-
import os

# Only include data dirs that actually exist on the build machine
datas = [('web', 'web')]

if os.path.isdir('data/reference'):
    datas.append(('data/reference', 'data/reference'))

if os.path.isdir('backend/models'):
    datas.append(('backend/models', 'backend/models'))

HIDDEN = [
    'uvicorn.protocols.http.h11_impl',
    'uvicorn.loop.asyncio',
    'uvicorn.lifespan.on',
    'uvicorn.lifespan.off',
    'uvicorn.loop.auto',
    'sklearn.ensemble._forest',
    'sklearn.utils._typedefs',
    'fastapi',
    'pydantic',
    'cobra',
    'depinfo',
    'tkinter',
    'tkinter.ttk',
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

a = Analysis(
    ['launcher.pyw'],          # 入口改为无窗口启动器
    pathex=['backend'],
    binaries=[],
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
