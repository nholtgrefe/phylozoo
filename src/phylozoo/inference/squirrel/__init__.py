"""
Squirrel submodule.

This submodule provides classes for working with quarnets and quarnetsets.
"""

from .quarnet import Quarnet, CycleQuarnet, QuartetTree, SplitQuarnet, SingleTriangle
from .quarnetset import QuarnetSet
from .tstar_tree import bstar, tstar_tree

__all__ = [
    "Quarnet",
    "CycleQuarnet",
    "QuartetTree",
    "SplitQuarnet",
    "SingleTriangle",
    "QuarnetSet",
    "bstar",
    "tstar_tree",
]

