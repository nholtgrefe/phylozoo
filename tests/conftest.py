"""
Pytest configuration and shared fixtures.
"""

import sys
from pathlib import Path

# Add src directory to Python path for imports
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

import pytest


@pytest.fixture
def sample_split():
    """
    Create a sample Split for testing.

    Returns
    -------
    Split
        A sample split object
    """
    from phylozoo.core.structure import Split

    return Split({1, 2}, {3, 4})


@pytest.fixture
def sample_split_system():
    """
    Create a sample SplitSystem for testing.

    Returns
    -------
    SplitSystem
        A sample split system object
    """
    from phylozoo.core.structure import Split, SplitSystem

    split1 = Split({1, 2}, {3, 4})
    split2 = Split({1, 3}, {2, 4})
    return SplitSystem([split1, split2])


