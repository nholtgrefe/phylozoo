"""
Network generator domain exceptions for PhyloZoo.
"""

from __future__ import annotations

from phylozoo.utils.exceptions.general import PhyloZooValueError


class PhyloZooGeneratorError(PhyloZooValueError):
    """
    Base exception for network generator-related errors.
    """

    pass


class PhyloZooGeneratorStructureError(PhyloZooGeneratorError):
    """
    Raised when network generator structure is invalid.
    """

    pass


class PhyloZooGeneratorDegreeError(PhyloZooGeneratorError):
    """
    Raised when network generator node degree constraints are violated.
    """

    pass
