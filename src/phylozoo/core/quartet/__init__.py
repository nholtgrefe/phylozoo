"""
Quartet module.

This module provides classes for representing and working with quartets, which are
unrooted trees on 4 taxa. A quartet can either be resolved (with a single non-trivial
split) or unresolved (a star tree). The public API (Quartet, QuartetProfile,
QuartetProfileSet) is re-exported here; the implementation is split across the base,
qprofile, qprofileset, and qdistance submodules.
"""

from .base import Quartet
from .qprofile import QuartetProfile
from .qprofileset import QuartetProfileSet

__all__ = [
    "Quartet",
    "QuartetProfile",
    "QuartetProfileSet",
]

