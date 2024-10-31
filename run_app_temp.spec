# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

a = Analysis(
    ['run.py'],  # The main script to run
    pathex=[],   # Additional search paths for imports
    binaries=[],  # Additional binaries to include
    datas=[
        # Add paths to specific data files from the environment
        ("./myenv/Lib/site-packages/altair/vegalite/v5/schema/vega-lite-schema.json", "./altair/vegalite/v5/schema/"),
        ("./myenv/Lib/site-packages/streamlit", "./streamlit/"),
        ("./myenv/Lib/site-packages/importlib_metadata", "./importlib_metadata/"),
        ("./myenv/Lib/site-packages/streamlit/static", "./streamlit/static"),
        ("./myenv/Lib/site-packages/streamlit/runtime", "./streamlit/runtime"),
        ("./myenv/Lib/site-packages/plotly", "./plotly/"),
        ("./myenv/Lib/site-packages/pingouin", "./pingouin/"),
        ("./myenv/Lib/site-packages/kaleido", "./kaleido/"),
        ("./myenv/Lib/site-packages/openpyxl", "./openpyxl/"),
        ("./myenv/Lib/site-packages/scikit_posthocs", "./scikit_posthocs/"),
        ("./myenv/Lib/site-packages/scikit_bio", "./scikit_bio/"),  # Added scikit-bio
        ("./myenv/Lib/site-packages/gnpsdata", "./gnpsdata/"),
        ("./myenv/Lib/site-packages/sklearn", "./sklearn/"),
        ("./myenv/Lib/site-packages/networkx", "./networkx/"),
        ("./myenv/Lib/site-packages/tabulate", "./tabulate/"),
        ("./myenv/Lib/site-packages/pandas_flavor", "./pandas_flavor/"),
        ("./myenv/Lib/site-packages/numpy", "./numpy/"),
        ("./myenv/Lib/site-packages/scipy", "./scipy/"),  # Added scipy
        ("./myenv/Lib/site-packages/streamlit_agraph", "./streamlit_agraph/"),
        ("./myenv/Lib/site-packages/pygraphviz", "./pygraphviz/"),
    ],
    hiddenimports=[],  # Specify hidden imports if needed
    hookspath=['./hooks'],  # Custom hook path if you have specific hooks
    hooksconfig={},
    runtime_hooks=[],  # Additional runtime hooks, if needed
    excludes=[],  # Packages to exclude from the build
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,  # Use encryption for the archive if needed
    noarchive=False,  # Store additional files in the archive
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='ChemProp2-App',  # Name of the executable
    debug=False,  # Set to True to enable debug logs
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,  # Use UPX compression
    upx_exclude=[],  # Exclude specific binaries from UPX compression if needed
    runtime_tmpdir=None,  # Temporary directory for runtime
    console=True,  # Set to True if it's a console app, False for a windowed app
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
