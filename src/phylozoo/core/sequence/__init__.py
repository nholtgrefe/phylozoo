"""
Sequence module.

This module provides classes and functions for working with multiple sequence alignments (MSAs). Sequences are stored internally as NumPy
arrays for efficient computation. The public API (MSA, bootstrap, bootstrap_per_gene,
hamming_distances) is re-exported here; the implementation is split across the
base, distances, bootstrap, and io submodules.
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

