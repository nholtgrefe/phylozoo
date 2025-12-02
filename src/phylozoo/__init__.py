"""
PhyloZoo: A phylogenetic analysis package.

This package provides tools for working with phylogenetic networks, trees,
and related structures.
"""

__version__ = "0.1.0"

# Import structures fundamental classes for convenient access
from .structures import (
    Split,
    SplitSet,
    SplitSystem,
    QuartetSplitSet,
    Partition,
    CircularOrdering,
    CircularSetOrdering,
)

# Import network classes
from .networks import (
    DirectedNetwork,
    SemiDirectedNetwork,
    random_semi_directed_network,
    MixedGraph,
    MultiMixedGraph,
)

__all__ = [
    "__version__",
    # Structures
    "Split",
    "SplitSet",
    "SplitSystem",
    "QuartetSplitSet",
    "Partition",
    "CircularOrdering",
    "CircularSetOrdering",
    # Networks
    "DirectedNetwork",
    "SemiDirectedNetwork",
    "random_semi_directed_network",
    "MixedGraph",
    "MultiMixedGraph",
]
