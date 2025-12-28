"""Pytest configuration for TSQVT tests."""

import sys
from pathlib import Path

# Add src to path for all tests
src_path = Path(__file__).parent.parent / 'src'
sys.path.insert(0, str(src_path))
