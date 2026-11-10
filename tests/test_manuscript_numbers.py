"""Every number quoted in the manuscript is recomputed from source_data/ and matches the text."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_manuscript_numbers_match():
    r = subprocess.run([sys.executable, str(ROOT / "src/02_statistical_tests/verify_manuscript_numbers.py")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:]
