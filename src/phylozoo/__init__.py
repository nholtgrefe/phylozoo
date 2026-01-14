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
    # Structures
    Split,
    SplitSystem,
    # Distance
    DistanceMatrix,
    # Sequence
    MSA,
    # Quartet
    Quartet,
    QuartetProfile,
    QuartetProfileSet,
)

# Import alias modules to register them
# These allow: from phylozoo.core.dnetwork import ... instead of from phylozoo.core.network.dnetwork import ...
from .core import dnetwork  # noqa: F401
from .core import sdnetwork  # noqa: F401

# Import exceptions for top-level access
from .utils.exceptions import (
    PhyloZooError,
    PhyloZooNotImplementedError,
    PhyloZooValueError,
    PhyloZooTypeError,
    PhyloZooRuntimeError,
    PhyloZooImportError,
    PhyloZooAttributeError,
    PhyloZooNetworkError,
    PhyloZooNetworkStructureError,
    PhyloZooNetworkDegreeError,
    PhyloZooNetworkAttributeError,
    PhyloZooIOError,
    PhyloZooParseError,
    PhyloZooFormatError,
    PhyloZooAlgorithmError,
    PhyloZooVisualizationError,
    PhyloZooLayoutError,
    PhyloZooBackendError,
    PhyloZooStateError,
    PhyloZooWarning,
    PhyloZooIdentifierWarning,
    PhyloZooEmptyNetworkWarning,
    PhyloZooSingleNodeNetworkWarning,
    PhyloZooGeneratorError,
    PhyloZooGeneratorStructureError,
    PhyloZooGeneratorDegreeError,
)

__all__ = [
    "__version__",
    # Core - Networks
    "DirectedPhyNetwork",
    "SemiDirectedPhyNetwork",
    # Core - Structures
    "Split",
    "SplitSystem",
    # Core - Distance
    "DistanceMatrix",
    # Core - Sequence
    "MSA",
    # Core - Quartet
    "Quartet",
    "QuartetProfile",
    "QuartetProfileSet",
    # Exceptions
    "PhyloZooError",
    "PhyloZooNotImplementedError",
    "PhyloZooValueError",
    "PhyloZooTypeError",
    "PhyloZooRuntimeError",
    "PhyloZooImportError",
    "PhyloZooAttributeError",
    "PhyloZooNetworkError",
    "PhyloZooNetworkStructureError",
    "PhyloZooNetworkDegreeError",
    "PhyloZooNetworkAttributeError",
    "PhyloZooIOError",
    "PhyloZooParseError",
    "PhyloZooFormatError",
    "PhyloZooAlgorithmError",
    "PhyloZooVisualizationError",
    "PhyloZooLayoutError",
    "PhyloZooBackendError",
    "PhyloZooStateError",
    "PhyloZooWarning",
    "PhyloZooIdentifierWarning",
    "PhyloZooEmptyNetworkWarning",
    "PhyloZooSingleNodeNetworkWarning",
    "PhyloZooGeneratorError",
    "PhyloZooGeneratorStructureError",
    "PhyloZooGeneratorDegreeError",
]
