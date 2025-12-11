"""
Validation utilities for classes with state validation.

This module provides decorators and context managers for validating object state
after mutations, with support for nested validation contexts and cache clearing.
"""

from contextlib import contextmanager
from functools import wraps
from typing import Any, Callable, TypeVar

T = TypeVar('T')


def validates(validate_method_name: str = '_validate') -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator factory to automatically validate object state after mutation methods.
    
    Automatically calls the validation method after the decorated method executes,
    but only if validation mode is enabled (respects no_validation context manager).
    
    Parameters
    ----------
    validate_method_name : str, optional
        Name of the validation method to call. By default '_validate'.
    
    Returns
    -------
    Callable
        Decorator that can be applied to mutation methods.
    
    Examples
    --------
    >>> class MyClass:
    ...     def __init__(self):
    ...         self._validate_mode = True
    ...         self._validation_depth = 0
    ...
    ...     @validates()
    ...     def add_item(self, item):
    ...         self.items.append(item)
    ...         # Automatically calls _validate() if validation enabled
    ...
    ...     def _validate(self):
    ...         if len(self.items) > 10:
    ...             raise ValueError("Too many items")
    
    >>> # Custom validation method name
    >>> class MyNetwork:
    ...     @validates('_validate_network')
    ...     def add_edge(self, u, v):
    ...         self.edges.append((u, v))
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(self: Any, *args: Any, **kwargs: Any) -> T:
            # Execute the mutation method
            result = func(self, *args, **kwargs)
            
            # Only validate if validation mode is enabled
            if getattr(self, '_validate_mode', True):
                if hasattr(self, validate_method_name):
                    validate_method = getattr(self, validate_method_name)
                    validate_method()
            
            return result
        
        return wrapper
    return decorator


@contextmanager
def no_validation_context(obj: Any, validate_method_name: str = '_validate'):
    """
    Context manager to disable validation for algorithms.
    
    Supports nested contexts - validation only occurs when exiting
    the outermost context.
    
    Parameters
    ----------
    obj : Any
        The object that supports validation.
    validate_method_name : str, optional
        Name of the validation method to call. By default '_validate'.
    
    Yields
    ------
    Any
        The object itself, for chaining operations.
    
    Examples
    --------
    >>> class MyClass:
    ...     def __init__(self):
    ...         self._validate_mode = True
    ...         self._validation_depth = 0
    ...
    ...     def no_validation(self):
    ...         return no_validation_context(self)
    ...
    >>> obj = MyClass()
    >>> with obj.no_validation():
    ...     obj.add_item(1)
    ...     with obj.no_validation():  # Nested
    ...         obj.add_item(2)
    ...     # Still in outer context - no validation yet
    ...     obj.add_item(3)
    ... # Now validates (exiting outermost context)
    """
    # Ensure validation state exists
    if not hasattr(obj, '_validation_depth'):
        obj._validation_depth = 0
    if not hasattr(obj, '_validate_mode'):
        obj._validate_mode = True
    
    obj._validation_depth += 1
    old_mode = obj._validate_mode
    obj._validate_mode = False
    
    try:
        yield obj
    finally:
        obj._validation_depth -= 1
        obj._validate_mode = old_mode
        
        # Only validate when exiting outermost context
        if obj._validation_depth == 0:
            if hasattr(obj, validate_method_name):
                validate_method = getattr(obj, validate_method_name)
                validate_method()

