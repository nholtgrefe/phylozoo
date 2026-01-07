"""
Directed multi-graph module.

This module provides the DirectedMultiGraph class.
"""

from .base import DirectedMultiGraph
from .isomorphism import is_isomorphic
from . import io

__all__ = [
    "DirectedMultiGraph",
    "is_isomorphic",
    "io",
]
