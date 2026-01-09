"""
Squirrel submodule.

This submodule provides classes for working with squirrel quartet profiles and profile sets.
"""

from .sqprofile import SqQuartetProfile
from .sqprofileset import SqQuartetProfileSet
from .qjoining import adapted_quartet_joining, quartet_joining
from .unresolve_tree import split_support, unresolve_tree
from .tstar_tree import bstar, tstar_tree
from .cycle_resolution import resolve_cycle_ordering

__all__ = [
    "SqQuartetProfile",
    "SqQuartetProfileSet",
    "adapted_quartet_joining",
    "quartet_joining",
    "split_support",
    "unresolve_tree",
    "bstar",
    "tstar_tree",
    "resolve_cycle_ordering",
]

