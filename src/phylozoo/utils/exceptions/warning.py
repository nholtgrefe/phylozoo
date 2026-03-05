"""
Warning classes for PhyloZoo.
"""

from __future__ import annotations


class PhyloZooWarning(UserWarning):
    """
    Base warning class for all PhyloZoo warnings.

    All custom warnings in PhyloZoo inherit from this class, allowing
    users to filter or catch all PhyloZoo-specific warnings. Inherits from
    UserWarning for compatibility with Python's warning system.

    Examples
    --------
    >>> import warnings
    >>> from phylozoo.utils.exceptions import PhyloZooWarning
    >>>
    >>> warnings.filterwarnings('ignore', category=PhyloZooWarning)
    >>> warnings.warn("This is a PhyloZoo warning", PhyloZooWarning)
    """

    pass


class PhyloZooIdentifierWarning(PhyloZooWarning):
    """
    Warning for identifier-related issues.

    This warning is used when Python keywords are used as identifiers
    or when None is used as a value, which may cause unexpected behavior.
    """

    pass


class PhyloZooEmptyNetworkWarning(PhyloZooWarning):
    """
    Warning for empty network operations.

    This warning is used when an operation is attempted on an empty network.
    """

    pass


class PhyloZooSingleNodeNetworkWarning(PhyloZooWarning):
    """
    Warning for single-node network operations.

    This warning is used when an operation is attempted on a single-node network.
    """

    pass
