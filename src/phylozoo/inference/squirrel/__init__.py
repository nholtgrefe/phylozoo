"""
Squirrel submodule.

This submodule provides classes for working with quarnets and quarnetsets.
"""

from .quarnet import Quarnet, CycleQuarnet, QuartetTree, SplitQuarnet, SingleTriangle
from .quarnetset import QuarnetSet
from .qjoining import adapted_quartet_joining, quartet_joining
from .unresolve_tree import split_support, unresolve_tree
from .tstar_tree import bstar, tstar_tree

__all__ = [
    "Quarnet",
    "CycleQuarnet",
    "QuartetTree",
    "SplitQuarnet",
    "SingleTriangle",
    "QuarnetSet",
    "adapted_quartet_joining",
    "quartet_joining",
    "split_support",
    "unresolve_tree",
    "bstar",
    "tstar_tree",
]

