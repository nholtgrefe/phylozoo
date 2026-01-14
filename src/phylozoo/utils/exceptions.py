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


class PhyloZooValidationError(PhyloZooError, ValueError):
    """
    Raised when validation fails.
    
    This exception is used for general validation failures that don't fit
    into more specific categories. Inherits from both PhyloZooError and
    ValueError for backward compatibility.
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooValidationError
    >>> 
    >>> if not is_valid(data):
    ...     raise PhyloZooValidationError("Data validation failed")
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


class PhyloZooNetworkTopologyError(PhyloZooNetworkError):
    """
    Raised when network topology is invalid.
    
    This exception is used for topology-related errors such as:
    - Invalid blob structure
    - Invalid network type (e.g., level-k violations)
    - Invalid generator structure
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNetworkTopologyError
    >>> 
    >>> if network.level > max_level:
    ...     raise PhyloZooNetworkTopologyError(
    ...         f"Network level {network.level} exceeds maximum {max_level}"
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


class PhyloZooFileNotFoundError(PhyloZooIOError, FileNotFoundError):
    """
    Raised when a file is not found.
    
    This exception is used when attempting to read or access a file that
    doesn't exist. Inherits from both PhyloZooIOError and FileNotFoundError
    for backward compatibility.
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooFileNotFoundError
    >>> 
    >>> if not path.exists():
    ...     raise PhyloZooFileNotFoundError(f"File not found: {path}")
    """
    pass


class PhyloZooParseError(PhyloZooError):
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


class PhyloZooInferenceError(PhyloZooAlgorithmError):
    """
    Raised when inference algorithms fail.
    
    This exception is used for inference algorithm failures such as:
    - Network reconstruction failures
    - Profile set inference errors
    - Tree inference errors
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooInferenceError
    >>> 
    >>> if not can_reconstruct(profileset):
    ...     raise PhyloZooInferenceError(
    ...         "Cannot reconstruct network from profile set"
    ...     )
    """
    pass


class PhyloZooComputationError(PhyloZooAlgorithmError):
    """
    Raised when computation fails.
    
    This exception is used for general computation failures such as:
    - Distance matrix computation errors
    - Diversity calculation errors
    - Optimization failures
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooComputationError
    >>> 
    >>> if computation_failed:
    ...     raise PhyloZooComputationError(
    ...         "Computation failed due to numerical instability"
    ...     )
    """
    pass


# ============================================================================
# State/Operation Domain Exceptions
# ============================================================================

class PhyloZooStateError(PhyloZooError):
    """
    Raised when an operation is attempted in an invalid state.
    
    This exception is used when an operation requires a certain state
    (e.g., object must be initialized, plot must be created) but that
    state is not met.
    
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


class PhyloZooNotFoundError(PhyloZooError):
    """
    Raised when an entity is not found.
    
    This exception is used when looking up entities that don't exist,
    such as taxa, nodes, profiles, or other identifiers.
    
    Examples
    --------
    >>> from phylozoo.utils.exceptions import PhyloZooNotFoundError
    >>> 
    >>> node_id = network.get_node_id(taxon)
    >>> if node_id is None:
    ...     raise PhyloZooNotFoundError(
    ...         f"Taxon '{taxon}' not found in network"
    ...     )
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
