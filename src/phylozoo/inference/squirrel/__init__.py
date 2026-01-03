"""
Squirrel submodule.

This submodule provides classes for working with quarnets and quarnetsets.
"""

from .quarnet import Quarnet, CycleQuarnet, QuartetTree, SplitQuarnet, SingleTriangle
from .quarnetset import QuarnetSet
from .bstar import bstar

__all__ = [
    "Quarnet",
    "CycleQuarnet",
    "QuartetTree",
    "SplitQuarnet",
    "SingleTriangle",
    "QuarnetSet",
    "bstar",
]

