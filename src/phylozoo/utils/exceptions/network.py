"""
Network domain exceptions for PhyloZoo.
"""

from __future__ import annotations

from phylozoo.utils.exceptions.general import PhyloZooValueError


class PhyloZooNetworkError(PhyloZooValueError):
    """
    Base exception for network-related errors.

    All network-specific errors inherit from this class, allowing users
    to catch all network errors with a single except clause.
    """

    pass


class PhyloZooNetworkStructureError(PhyloZooNetworkError):
    """
    Raised when network structure is invalid.

    This exception is used for structural validation failures such as:
    - Directed cycles
    - Disconnected networks
    - Self-loops
    - Invalid connectivity
    """

    pass


class PhyloZooNetworkDegreeError(PhyloZooNetworkError):
    """
    Raised when node degree constraints are violated.

    This exception is used when nodes have invalid degrees, such as:
    - Leaf nodes with wrong in-degree or out-degree
    - Internal nodes with invalid degree combinations
    - Root nodes with wrong in-degree
    """

    pass


class PhyloZooNetworkAttributeError(PhyloZooNetworkError):
    """
    Raised when network edge attributes are invalid.

    This exception is used for attribute validation failures such as:
    - Invalid gamma values (not in [0, 1] or don't sum to 1.0)
    - Invalid bootstrap values (not in [0, 1])
    - Invalid branch_length values (inconsistent across parallel edges)
    """

    pass
