"""
Algorithm domain exceptions for PhyloZoo.
"""

from __future__ import annotations

from phylozoo.utils.exceptions.base import PhyloZooError


class PhyloZooAlgorithmError(PhyloZooError):
    """
    Base exception for algorithm-related errors.

    All algorithm-specific errors inherit from this class, allowing users
    to catch all algorithm errors with a single except clause.
    """

    pass
