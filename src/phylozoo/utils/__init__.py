"""
Utility modules for PhyloZoo.

This package contains various utility functions and classes.
"""

from .exceptions import (
    PhyloZooError,
    PhyloZooNotImplementedError,
    PhyloZooValueError,
    PhyloZooTypeError,
    PhyloZooNetworkError,
    PhyloZooNetworkStructureError,
    PhyloZooNetworkDegreeError,
    PhyloZooNetworkAttributeError,
    PhyloZooIOError,
    PhyloZooImportError,
    PhyloZooParseError,
    PhyloZooFormatError,
    PhyloZooAlgorithmError,
    PhyloZooStateError,
    PhyloZooVisualizationError,
    PhyloZooLayoutError,
    PhyloZooBackendError,
    PhyloZooWarning,
    PhyloZooIdentifierWarning,
)
from .identifier_warnings import warn_on_keyword, warn_on_none_value
from .validation import no_validation, validation_aware

__all__ = [
    # Exceptions
    "PhyloZooError",
    "PhyloZooNotImplementedError",
    "PhyloZooValueError",
    "PhyloZooTypeError",
    "PhyloZooNetworkError",
    "PhyloZooNetworkStructureError",
    "PhyloZooNetworkDegreeError",
    "PhyloZooNetworkAttributeError",
    "PhyloZooIOError",
    "PhyloZooImportError",
    "PhyloZooParseError",
    "PhyloZooFormatError",
    "PhyloZooAlgorithmError",
    "PhyloZooStateError",
    "PhyloZooVisualizationError",
    "PhyloZooLayoutError",
    "PhyloZooBackendError",
    # Warnings
    "PhyloZooWarning",
    "PhyloZooIdentifierWarning",
    # Other utilities
    "warn_on_keyword",
    "warn_on_none_value",
    "no_validation",
    "validation_aware",
]
