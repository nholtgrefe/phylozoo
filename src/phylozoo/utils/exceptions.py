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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNotImplementedError
    >>> 
    >>> def some_method(self):
    ...     raise PhyloZooNotImplementedError(
    ...         "This method is not yet implemented"
    ...     )
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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooValueError
    >>> 
    >>> # Value out of range
    >>> if k < 0:
    ...     raise PhyloZooValueError(f"k must be non-negative, got {k}")
    >>> 
    >>> # Invalid value (right type, wrong value)
    >>> if bootstrap < 0.0 or bootstrap > 1.0:
    ...     raise PhyloZooValueError(f"bootstrap must be in [0.0, 1.0], got {bootstrap}")
    """
    pass


class PhyloZooTypeError(PhyloZooError, TypeError):
    """
    Raised when a type is incorrect.
    
    This exception is used when an argument has the wrong type (e.g., string
    when number expected, wrong class instance). Inherits from both PhyloZooError
    and TypeError for backward compatibility.
    
    Note: For value errors (right type, wrong value), use PhyloZooValueError instead.
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooTypeError
    >>> 
    >>> # Wrong type
    >>> if not isinstance(obj, ExpectedType):
    ...     raise PhyloZooTypeError(
    ...         f"Expected {ExpectedType}, got {type(obj)}"
    ...     )
    >>> 
    >>> # Type check (e.g., "must be numeric")
    >>> if not isinstance(value, (int, float)):
    ...     raise PhyloZooTypeError(
    ...         f"Expected numeric type, got {type(value).__name__}"
    ...     )
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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooImportError
    >>> 
    >>> try:
    ...     import optional_module
    ... except ImportError as e:
    ...     raise PhyloZooImportError(
    ...         f"Optional module 'optional_module' is required. "
    ...         "Install with: pip install optional_module"
    ...     ) from e
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

class PhyloZooNetworkError(PhyloZooError):
    """
    Base exception for network-related errors.
    
    All network-specific errors inherit from this class, allowing users
    to catch all network errors with a single except clause.
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNetworkError
    >>> 
    >>> try:
    ...     network.validate()
    ... except PhyloZooNetworkError as e:
    ...     print(f"Network error: {e}")
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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNetworkStructureError
    >>> 
    >>> if has_cycles(network):
    ...     raise PhyloZooNetworkStructureError(
    ...         "Network contains directed cycles"
    ...     )
    """
    pass


class PhyloZooNetworkDegreeError(PhyloZooNetworkError):
    """
    Raised when node degree constraints are violated.
    
    This exception is used when nodes have invalid degrees, such as:
    - Leaf nodes with wrong in-degree or out-degree
    - Internal nodes with invalid degree combinations
    - Root nodes with wrong in-degree
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNetworkDegreeError
    >>> 
    >>> if leaf.indegree() != 1:
    ...     raise PhyloZooNetworkDegreeError(
    ...         f"Leaf node {leaf} has invalid in-degree"
    ...     )
    """
    pass


class PhyloZooNetworkAttributeError(PhyloZooNetworkError):
    """
    Raised when network edge attributes are invalid.
    
    This exception is used for attribute validation failures such as:
    - Invalid gamma values (not in [0, 1] or don't sum to 1.0)
    - Invalid bootstrap values (not in [0, 1])
    - Invalid branch_length values (inconsistent across parallel edges)
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNetworkAttributeError
    >>> 
    >>> if gamma_sum != 1.0:
    ...     raise PhyloZooNetworkAttributeError(
    ...         f"Gamma values must sum to 1.0, got {gamma_sum}"
    ...     )
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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooIOError
    >>> 
    >>> try:
    ...     content = file.read()
    ... except IOError as e:
    ...     raise PhyloZooIOError(f"Failed to read file: {e}") from e
    """
    pass


class PhyloZooParseError(PhyloZooIOError):
    """
    Raised when parsing fails.
    
    This exception is used for parsing errors in various formats such as:
    - eNewick parsing errors
    - Newick parsing errors
    - Other format parsing errors
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooParseError
    >>> 
    >>> try:
    ...     parse_enewick(string)
    ... except Exception as e:
    ...     raise PhyloZooParseError(f"Failed to parse eNewick: {e}") from e
    """
    pass


class PhyloZooFormatError(PhyloZooIOError):
    """
    Raised when format-related operations fail.
    
    This exception is used for format-related errors such as:
    - Unsupported formats
    - Format conversion failures
    - Format detection failures
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooFormatError
    >>> 
    >>> if format not in supported_formats:
    ...     raise PhyloZooFormatError(
    ...         f"Format '{format}' not supported"
    ...     )
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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooAlgorithmError
    >>> 
    >>> try:
    ...     result = algorithm.compute()
    ... except PhyloZooAlgorithmError as e:
    ...     print(f"Algorithm error: {e}")
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
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooVisualizationError
    >>> 
    >>> try:
    ...     plot_network(network)
    ... except PhyloZooVisualizationError as e:
    ...     print(f"Visualization error: {e}")
    """
    pass


class PhyloZooLayoutError(PhyloZooVisualizationError):
    """
    Raised when layout computation fails.
    
    This exception is used for layout computation errors such as:
    - Empty network/graph layout errors
    - Invalid layout algorithm
    - Layout computation failures
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooLayoutError
    >>> 
    >>> if network.number_of_nodes() == 0:
    ...     raise PhyloZooLayoutError(
    ...         "Cannot compute layout for empty network"
    ...     )
    """
    pass


class PhyloZooBackendError(PhyloZooVisualizationError):
    """
    Raised when backend operations fail.
    
    This exception is used for backend-related errors such as:
    - Backend not registered
    - Backend initialization failures
    - Backend operation failures
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooBackendError
    >>> 
    >>> if backend_name not in registered_backends:
    ...     raise PhyloZooBackendError(
    ...         f"Backend '{backend_name}' not registered"
    ...     )
    """
    pass


class PhyloZooStateError(PhyloZooVisualizationError):
    """
    Raised when an operation is attempted in an invalid state.
    
    This exception is used when an operation requires a certain state
    (e.g., plot must be created, figure must be initialized) but that
    state is not met. This is specific to visualization operations.
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooStateError
    >>> 
    >>> if self._figure is None:
    ...     raise PhyloZooStateError(
    ...         "Figure must be created before showing"
    ...     )
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
    
    Examples
    --------
    >>> import warnings
    >>> from phylozoo.utils.exceptions import PhyloZooIdentifierWarning
    >>> 
    >>> warnings.warn(
    ...     "Identifier 'for' is a Python keyword",
    ...     PhyloZooIdentifierWarning
    ... )
    """
    pass
