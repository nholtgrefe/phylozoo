"""
General-purpose exceptions for PhyloZoo.
"""

from __future__ import annotations

from phylozoo.utils.exceptions.base import PhyloZooError


class PhyloZooNotImplementedError(PhyloZooError, NotImplementedError):
    """
    Raised when a feature is not implemented.

    This exception is used when a method or feature is declared but not
    yet implemented. Inherits from both PhyloZooError and NotImplementedError
    for backward compatibility.
    """

    pass


class PhyloZooValueError(PhyloZooError, ValueError):
    """
    Raised when a value is invalid or out of range.

    This exception is used for parameter validation errors where the type is
    correct but the value is inappropriate (e.g., out of range, negative when
    positive expected). Inherits from both PhyloZooError and ValueError for
    backward compatibility.

    Note: For type errors (wrong type), use PhyloZooTypeError instead.
    """

    pass


class PhyloZooTypeError(PhyloZooError, TypeError):
    """
    Raised when a type is incorrect.

    This exception is used when an argument has the wrong type (e.g., string
    when number expected, wrong class instance). Inherits from both PhyloZooError
    and TypeError for backward compatibility.

    Note: For value errors (right type, wrong value), use PhyloZooValueError instead.
    """

    pass


class PhyloZooRuntimeError(PhyloZooError, RuntimeError):
    """
    Raised when a runtime error occurs.

    This exception is used for general runtime errors. Inherits from both
    PhyloZooError and RuntimeError for backward compatibility.
    """

    pass


class PhyloZooImportError(PhyloZooError, ImportError):
    """
    Raised when an import fails.

    This exception is used when a required module or dependency cannot be
    imported. Inherits from both PhyloZooError and ImportError for backward
    compatibility.
    """

    pass


class PhyloZooAttributeError(PhyloZooError, AttributeError):
    """
    Raised when an attribute is invalid.

    This exception is used when an attribute is invalid. Inherits from both
    PhyloZooError and AttributeError for backward compatibility.
    """

    pass
