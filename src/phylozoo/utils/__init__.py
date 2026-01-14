"""
Utility modules for PhyloZoo.

This package contains various utility functions and classes.
"""

from .exceptions import (
    PhyloZooError,
    PhyloZooNotImplementedError,
    PhyloZooValueError,
    PhyloZooTypeError,
    PhyloZooValidationError,
    PhyloZooNetworkError,
    PhyloZooNetworkStructureError,
    PhyloZooNetworkDegreeError,
    PhyloZooNetworkTopologyError,
    PhyloZooNetworkAttributeError,
    PhyloZooIOError,
    PhyloZooFileNotFoundError,
    PhyloZooImportError,
    PhyloZooParseError,
    PhyloZooFormatError,
    PhyloZooAlgorithmError,
    PhyloZooInferenceError,
    PhyloZooComputationError,
    PhyloZooStateError,
    PhyloZooNotFoundError,
    PhyloZooVisualizationError,
    PhyloZooLayoutError,
    PhyloZooBackendError,
)
from .identifier_warnings import warn_on_keyword, warn_on_none_value
from .validation import no_validation, validation_aware

__all__ = [
    # Exceptions
    "PhyloZooError",
    "PhyloZooNotImplementedError",
    "PhyloZooValueError",
    "PhyloZooTypeError",
    "PhyloZooValidationError",
    "PhyloZooNetworkError",
    "PhyloZooNetworkStructureError",
    "PhyloZooNetworkDegreeError",
    "PhyloZooNetworkTopologyError",
    "PhyloZooNetworkAttributeError",
    "PhyloZooIOError",
    "PhyloZooFileNotFoundError",
    "PhyloZooImportError",
    "PhyloZooParseError",
    "PhyloZooFormatError",
    "PhyloZooAlgorithmError",
    "PhyloZooInferenceError",
    "PhyloZooComputationError",
    "PhyloZooStateError",
    "PhyloZooNotFoundError",
    "PhyloZooVisualizationError",
    "PhyloZooLayoutError",
    "PhyloZooBackendError",
    # Other utilities
    "warn_on_keyword",
    "warn_on_none_value",
    "no_validation",
    "validation_aware",
]
