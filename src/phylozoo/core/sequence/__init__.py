"""
Sequence module.

This module provides classes for working with biological sequences.
"""

from .base import MSA
from .distances import hamming_distances

# Import I/O handlers to register them
from . import io  # noqa: F401

__all__ = [
    "MSA",
    "hamming_distances",
]

