"""
I/O domain exceptions for PhyloZoo.
"""

from __future__ import annotations

from phylozoo.utils.exceptions.base import PhyloZooError


class PhyloZooIOError(PhyloZooError, IOError):
    """
    Raised when I/O operations fail.

    This exception is used for general I/O errors. Inherits from both
    PhyloZooError and IOError=OSError for backward compatibility.
    """

    pass


class PhyloZooParseError(PhyloZooIOError):
    """
    Raised when parsing fails.

    This exception is used for parsing errors in various formats such as:

    - eNewick parsing errors
    - Newick parsing errors
    - Other format parsing errors
    """

    pass


class PhyloZooFormatError(PhyloZooIOError):
    """
    Raised when format-related operations fail.

    This exception is used for format-related errors such as:

    - Unsupported formats
    - Format conversion failures
    - Format detection failures
    """

    pass
