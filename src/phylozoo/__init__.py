"""
PhyloZoo: A phylogenetic analysis package.

This package provides tools for working with phylogenetic networks, trees,
and related structures.
"""

__version__ = "0.1.0"

# Import core classes for convenient access
from .core import (
    # Networks
    DirectedPhyNetwork,
    SemiDirectedPhyNetwork,
    MixedMultiGraph,
    DirectedMultiGraph,
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

# Make core submodules available at package level
from . import core
# Import submodules so they're accessible as phylozoo.network, phylozoo.primitives, etc.
from .core import network as _network_module
from .core import primitives as _primitives_module
from .core import split as _split_module
from .core import distance as _distance_module
from .core import sequence as _sequence_module

# Re-export submodules at package level
network = _network_module
primitives = _primitives_module
split = _split_module
distance = _distance_module
sequence = _sequence_module

# Import inference classes
from .inference import (
    NetworkInferrer,
    infer_network_from_msa,
)

__all__ = [
    "__version__",
    # Core module
    "core",
    # Core submodules (available as phylozoo.network, phylozoo.primitives, etc.)
    "network",
    "primitives",
    "split",
    "distance",
    "sequence",
    # Core - Networks
    "DirectedPhyNetwork",
    "SemiDirectedPhyNetwork",
    "MixedMultiGraph",
    "DirectedMultiGraph",
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
