"""
PhyloZoo: A phylogenetic analysis package.

This package provides tools for working with phylogenetic networks, trees,
and related structures.
"""

__version__ = "0.1.0"

# Import core classes for convenient access
from .core import (
    # Networks
    DirectedNetwork,
    SemiDirectedNetwork,
    random_semi_directed_network,
    MixedGraph,
    MultiMixedGraph,
    # Structures
    Split,
    SplitSystem,
    Partition,
    CircularOrdering,
    CircularSetOrdering,
    # Distance
    DistanceMatrix,
    # Sequence
    MSA,
)

# Import inference classes
from .inference import (
    NetworkInferrer,
    infer_network_from_msa,
)

__all__ = [
    "__version__",
    # Core - Networks
    "DirectedNetwork",
    "SemiDirectedNetwork",
    "random_semi_directed_network",
    "MixedGraph",
    "MultiMixedGraph",
    # Core - Structures
    "Split",
    "SplitSystem",
    "Partition",
    "CircularOrdering",
    "CircularSetOrdering",
    # Core - Distance
    "DistanceMatrix",
    # Core - Sequence
    "MSA",
    # Inference
    "NetworkInferrer",
    "infer_network_from_msa",
]
