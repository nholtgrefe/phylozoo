"""
Base exception for PhyloZoo.

All custom exceptions in PhyloZoo inherit from PhyloZooError.
"""

from __future__ import annotations


class PhyloZooError(Exception):
    """
    Base exception for all PhyloZoo errors.

    All custom exceptions in PhyloZoo inherit from this class, allowing
    users to catch all PhyloZoo-specific errors with a single except clause.

    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooError
    >>>
    >>> try:
    ...     # Some PhyloZoo operation
    ...     pass
    ... except PhyloZooError as e:
    ...     print(f"Caught PhyloZoo error: {e}")
    """

    pass
