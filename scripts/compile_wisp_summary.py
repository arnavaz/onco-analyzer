#!/usr/bin/env python3
"""
compile_wisp_summary.py
=======================
DEPRECATED: Replaced by compile_wgts_summary.py.
Forwards all arguments to compile_wgts_summary.py for backward compatibility.
"""
import sys
import subprocess
from pathlib import Path

script_dir = Path(__file__).resolve().parent
target_script = script_dir / "compile_wgts_summary.py"

if not target_script.is_file():
    print(f"Error: {target_script} not found.", file=sys.stderr)
    sys.exit(1)

cmd = [sys.executable, str(target_script)] + sys.argv[1:]
res = subprocess.run(cmd)
sys.exit(res.returncode)
