from pathlib import Path
import sys
import importlib.util

SCRIPT_DIR = Path(__file__).resolve().parent
TARGET = SCRIPT_DIR / "tools" / "sample_discovery.py"
sys.path.insert(0, str(TARGET.parent))

spec = importlib.util.spec_from_file_location("_sample_discovery_impl", TARGET)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
globals().update({name: value for name, value in module.__dict__.items() if not name.startswith("_")})


if __name__ == "__main__":
    module.main()
