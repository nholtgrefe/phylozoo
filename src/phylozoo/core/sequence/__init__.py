"""
Sequence module.

This module provides classes for working with biological sequences.
"""

from .base import MSA
from .bootstrap import bootstrap, bootstrap_per_gene
from .distances import hamming_distances

# Import I/O handlers to register them
from . import io  # noqa: F401

__all__ = [
    "MSA",
    "bootstrap",
    "bootstrap_per_gene",
    "hamming_distances",
]

