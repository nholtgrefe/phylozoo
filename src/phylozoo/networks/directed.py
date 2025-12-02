"""
Directed network module.

This module provides classes and functions for working with directed phylogenetic networks.
"""

from contextlib import contextmanager
from functools import wraps
from typing import Callable, List, Optional, Set, Tuple, TypeVar, Any

T = TypeVar('T')


def validates(func: Callable[..., T]) -> Callable[..., T]:
    """
    Decorator to automatically validate network state after mutation methods.
    
    Automatically calls _validate_network() after the method executes, but only
    if validation mode is enabled (respects no_validation context manager).
    
    This decorator should be applied to methods that mutate the network state.
    Methods that only return information should not use this decorator.
    
    Examples
    --------
    >>> class MyNetwork:
    ...     @validates
    ...     def add_edge(self, source, target):
    ...         self.edges.append((source, target))
    ...         # Automatically calls _validate_network() if validation enabled
    """
    @wraps(func)
    def wrapper(self: Any, *args: Any, **kwargs: Any) -> T:
        # Execute the mutation method
        result = func(self, *args, **kwargs)
        
        # Only validate if validation mode is enabled
        if getattr(self, '_validate_mode', True):
            if hasattr(self, '_validate_network'):
                self._validate_network()
        
        return result
    
    return wrapper


class DirectedNetwork:
    """
    A directed phylogenetic network.
    
    Supports nested validation contexts for efficient batch operations.
    """

    def __init__(
        self, nodes: Optional[Set[str]] = None, edges: Optional[List[Tuple[str, str]]] = None
    ) -> None:
        """
        Initialize a directed network.

        Parameters
        ----------
        nodes : Optional[Set[str]], optional
            Set of node identifiers, by default None
        edges : Optional[List[Tuple[str, str]]], optional
            List of directed edges as tuples, by default None
        """
        self.nodes: Set[str] = nodes or set()
        self.edges: List[Tuple[str, str]] = edges or []
        self._validate_mode: bool = True  # Default: validate
        self._validation_depth: int = 0  # Track nesting depth

    @contextmanager
    def no_validation(self):
        """
        Context manager to disable validation for algorithms.
        
        Supports nested contexts - validation only occurs when exiting
        the outermost context.
        
        Yields
        ------
        DirectedNetwork
            Self, for chaining operations
        
        Examples
        --------
        >>> network = DirectedNetwork()
        >>> with network.no_validation():
        ...     network.add_edge("A", "B")
        ...     with network.no_validation():  # Nested
        ...         network.add_edge("B", "C")
        ...     # Still in outer context - no validation yet
        ...     network.add_edge("C", "D")
        ... # Now validates (exiting outermost context)
        """
        self._validation_depth += 1
        old_mode = self._validate_mode
        self._validate_mode = False
        
        try:
            yield self
        finally:
            self._validation_depth -= 1
            self._validate_mode = old_mode
            
            # Only validate when exiting outermost context
            if self._validation_depth == 0:
                self._validate_network()

    @validates
    def add_node(self, node: str) -> None:
        """
        Add a node to the network.

        Parameters
        ----------
        node : str
            Node identifier to add
        """
        self.nodes.add(node)

    @validates
    def add_edge(self, source: str, target: str) -> None:
        """
        Add a directed edge to the network.

        Parameters
        ----------
        source : str
            Source node identifier
        target : str
            Target node identifier
        """
        if source not in self.nodes:
            self.add_node(source)
        if target not in self.nodes:
            self.add_node(target)
        self.edges.append((source, target))


    def _validate_network(self) -> None:
        """
        Validate entire network state.
        
        Checks that the network is a valid directed network (e.g., all edges
        reference existing nodes, network invariants are maintained, etc.).
        
        Raises
        ------
        ValueError
            If network is invalid
        """
        errors = []
        
        # Check all edges reference existing nodes
        for source, target in self.edges:
            if source not in self.nodes:
                errors.append(f"Edge source '{source}' not in nodes")
            if target not in self.nodes:
                errors.append(f"Edge target '{target}' not in nodes")
        
        # Add other network-wide validations as needed
        # - Check for self-loops if not allowed
        # - Check for cycles if not allowed
        # - Check node degrees
        # - Check network topology invariants
        # - etc.
        
        if errors:
            raise ValueError(f"Invalid network:\n" + "\n".join(errors))

    def __repr__(self) -> str:
        """
        Return string representation of the network.

        Returns
        -------
        str
            String representation
        """
        return f"DirectedNetwork(nodes={len(self.nodes)}, edges={len(self.edges)})"
