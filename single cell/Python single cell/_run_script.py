from __future__ import annotations

import runpy
import sys
from pathlib import Path


def run_relative(relative_path: str) -> None:
    script_dir = Path(__file__).resolve().parent
    target = script_dir / relative_path
    if not target.exists():
        raise FileNotFoundError("Target script not found: {}".format(target))

    sys.path.insert(0, str(target.parent))
    sys.path.insert(0, str(script_dir / "tools"))
    sys.path.insert(0, str(script_dir / "main"))
    sys.path.insert(0, str(script_dir / "subcluster"))
    sys.argv[0] = str(target)
    runpy.run_path(str(target), run_name="__main__")
