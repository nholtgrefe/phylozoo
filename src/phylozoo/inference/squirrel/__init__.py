"""
Squirrel submodule.

This submodule provides classes for working with squirrel quartet profiles and profile sets.
"""

from .sqprofile import SqQuartetProfile
from .sqprofileset import SqQuartetProfileSet
from .qjoining import adapted_quartet_joining, quartet_joining
from .unresolve_tree import split_support, unresolve_tree
from .tstar_tree import bstar, tstar_tree
from .cycle_resolution import (
    _qprofiles_to_circular_ordering,
    _qprofiles_to_hybrid_ranking,
    _insert_cycle,
    resolve_cycles,
)
from .qsimilarity import sqprofileset_similarity, sqprofileset_from_network
from .delta_heuristic import delta_heuristic
from .squirrel import squirrel

__all__ = [
    "SqQuartetProfile",
    "SqQuartetProfileSet",
    "adapted_quartet_joining",
    "quartet_joining",
    "split_support",
    "unresolve_tree",
    "bstar",
    "tstar_tree",
    "_qprofiles_to_circular_ordering",
    "_qprofiles_to_hybrid_ranking",
    "_insert_cycle",
    "resolve_cycles",
    "sqprofileset_similarity",
    "sqprofileset_from_network",
    "delta_heuristic",
    "squirrel",
]

