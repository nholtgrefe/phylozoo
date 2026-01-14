"""
Custom exception hierarchy for PhyloZoo.

This module provides a comprehensive set of custom exceptions organized by domain.
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


# ============================================================================
# Base Exception
# ============================================================================

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


# ============================================================================
# General Purpose Exceptions
# ============================================================================

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

# ============================================================================
# Network Domain Exceptions
# ============================================================================

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


# ============================================================================
# I/O Domain Exceptions
# ============================================================================

class PhyloZooIOError(PhyloZooError, IOError):
    """
    Raised when I/O operations fail.
    
    This exception is used for general I/O errors. Inherits from both
    PhyloZooError and IOError for backward compatibility.
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


# ============================================================================
# Algorithm Domain Exceptions
# ============================================================================

class PhyloZooAlgorithmError(PhyloZooError):
    """
    Base exception for algorithm-related errors.
    
    All algorithm-specific errors inherit from this class, allowing users
    to catch all algorithm errors with a single except clause.
    """
    pass


# ============================================================================
# Visualization Domain Exceptions
# ============================================================================

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


# ============================================================================
# Warning Classes
# ============================================================================

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

###############################################################################
# Network Generator Domain Exceptions
###############################################################################

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