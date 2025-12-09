"""
Pytest configuration and shared fixtures for sdnetwork tests.
"""

from contextlib import contextmanager

import pytest


@contextmanager
def expect_mixed_network_warning():
    """
    Context manager to catch the expected MixedPhyNetwork validity warning.
    
    This warning is raised during network initialization and when validate()
    is called explicitly.
    """
    with pytest.warns(UserWarning, match="Validity is not fully checked for MixedPhyNetworks"):
        yield



