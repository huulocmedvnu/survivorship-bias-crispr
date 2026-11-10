"""Repository paths. Every script resolves files from the repository root, so they run from any directory."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"                 # raw public data (not distributed, see README)
SOURCE = ROOT / "source_data"        # processed per-perturbation tables (distributed)
FIGURES = ROOT / "output" / "figures"
TABLES = ROOT / "output" / "tables"
SCORER = ROOT / "rust_scorer" / "target" / "release" / "survb"
