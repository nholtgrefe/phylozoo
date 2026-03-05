"""
Visualization domain exceptions for PhyloZoo.
"""

from __future__ import annotations

from phylozoo.utils.exceptions.base import PhyloZooError


class PhyloZooVisualizationError(PhyloZooError):
    """
    Base exception for visualization-related errors.

    All visualization-specific errors inherit from this class, allowing
    users to catch all visualization errors with a single except clause.
    """

    pass


class PhyloZooLayoutError(PhyloZooVisualizationError):
    """
    Raised when layout computation fails.

    This exception is used for layout computation errors such as:
    - Empty network/graph layout errors
    - Invalid layout algorithm
    - Layout computation failures
    """

    pass


class PhyloZooBackendError(PhyloZooVisualizationError):
    """
    Raised when backend operations fail.

    This exception is used for backend-related errors such as:
    - Backend not registered
    - Backend initialization failures
    - Backend operation failures
    """

    pass


class PhyloZooStateError(PhyloZooVisualizationError):
    """
    Raised when an operation is attempted in an invalid state.

    This exception is used when an operation requires a certain state
    (e.g., plot must be created, figure must be initialized) but that
    state is not met. This is specific to visualization operations.
    """

    pass
