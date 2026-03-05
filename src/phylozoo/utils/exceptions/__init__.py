"""
Custom exception hierarchy for PhyloZoo.

This package provides a comprehensive set of custom exceptions organized by domain.
All exceptions inherit from PhyloZooError, allowing easy catching of PhyloZoo-specific
errors. Many also inherit from built-in exceptions (e.g., ValueError) for backward
compatibility.

Examples
--------
>>> from phylozoo.utils.exceptions import PhyloZooError, PhyloZooNetworkError
>>>
>>> # Catch all PhyloZoo errors
>>> try:
...     network.validate()
... except PhyloZooError as e:
...     print(f"PhyloZoo error: {e}")
>>>
>>> # Catch specific domain errors
>>> try:
...     network.validate()
... except PhyloZooNetworkError as e:
...     print(f"Network error: {e}")
>>>
>>> # Still works with built-in exceptions (backward compatible)
>>> try:
...     network.validate()
... except ValueError as e:
...     print(f"Value error: {e}")
"""

from __future__ import annotations

from phylozoo.utils.exceptions.algorithm import PhyloZooAlgorithmError
from phylozoo.utils.exceptions.base import PhyloZooError
from phylozoo.utils.exceptions.general import (
    PhyloZooAttributeError,
    PhyloZooImportError,
    PhyloZooNotImplementedError,
    PhyloZooRuntimeError,
    PhyloZooTypeError,
    PhyloZooValueError,
)
from phylozoo.utils.exceptions.generator import (
    PhyloZooGeneratorDegreeError,
    PhyloZooGeneratorError,
    PhyloZooGeneratorStructureError,
)
from phylozoo.utils.exceptions.io import (
    PhyloZooFormatError,
    PhyloZooIOError,
    PhyloZooParseError,
)
from phylozoo.utils.exceptions.network import (
    PhyloZooNetworkAttributeError,
    PhyloZooNetworkDegreeError,
    PhyloZooNetworkError,
    PhyloZooNetworkStructureError,
)
from phylozoo.utils.exceptions.utils import warn_on_keyword, warn_on_none_value
from phylozoo.utils.exceptions.visualization import (
    PhyloZooBackendError,
    PhyloZooLayoutError,
    PhyloZooStateError,
    PhyloZooVisualizationError,
)
from phylozoo.utils.exceptions.warning import (
    PhyloZooEmptyNetworkWarning,
    PhyloZooIdentifierWarning,
    PhyloZooSingleNodeNetworkWarning,
    PhyloZooWarning,
)

__all__ = [
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
    "warn_on_keyword",
    "warn_on_none_value",
]
