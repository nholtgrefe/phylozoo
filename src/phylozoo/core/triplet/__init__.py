"""
Triplet module.

This module provides classes for representing and working with triplets, which are
rooted trees on 3 taxa. A triplet can either be resolved (with a single trivial 1|2
split a|bc, where ``a`` is a child of the root and ``b, c`` form a cherry under an
internal node) or unresolved (a star tree where the root has three direct leaf
children). The public API (Triplet, TripletProfile, TripletProfileSet) is re-exported
here; the implementation is split across the base, tprofile, and tprofileset
submodules.
"""

from .base import Triplet
from .tprofile import TripletProfile
from .tprofileset import TripletProfileSet

__all__ = [
    "Triplet",
    "TripletProfile",
    "TripletProfileSet",
]
