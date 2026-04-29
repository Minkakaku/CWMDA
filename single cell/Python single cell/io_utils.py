from pathlib import Path
import sys
import importlib.util

SCRIPT_DIR = Path(__file__).resolve().parent
TARGET = SCRIPT_DIR / "tools" / "io_utils.py"
sys.path.insert(0, str(TARGET.parent))

spec = importlib.util.spec_from_file_location("_io_utils_impl", TARGET)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
globals().update({name: value for name, value in module.__dict__.items() if not name.startswith("_")})
