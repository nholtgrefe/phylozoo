"""
Caching utilities for computed properties.

This module provides decorators for caching computed properties that depend on
mutable object state. The cache is automatically cleared when the object is
modified through designated mutation methods.
"""

from functools import wraps
from typing import Any, Callable, TypeVar

T = TypeVar('T')


def cached_property(func: Callable[[Any], T]) -> property:
    """
    Cached property decorator that can be cleared when object state changes.
    
    This decorator caches the result of a property method. The cache is
    automatically cleared when mutation methods decorated with `clears_cache`
    are called.
    
    Parameters
    ----------
    func : Callable[[Any], T]
        The property method to cache.
    
    Returns
    -------
    property
        A cached property that automatically clears when object state changes.
    
    Examples
    --------
    >>> class MyClass:
    ...     def __init__(self):
    ...         self._data = [1, 2, 3]
    ...         self._cache_version = 0
    ...
    ...     @cached_property
    ...     def expensive_computation(self):
    ...         return sum(self._data) ** 2
    ...
    ...     @clears_cache
    ...     def add_item(self, item):
    ...         self._data.append(item)
    """
    cache_name = f'_cache_{func.__name__}'
    version_name = f'_cache_version_{func.__name__}'
    
    @property
    @wraps(func)
    def wrapper(self: Any) -> T:
        """
        Wrapper that checks cache validity and computes if needed.
        
        Parameters
        ----------
        self : Any
            The instance containing the property.
        
        Returns
        -------
        T
            The cached or newly computed value.
        """
        # Ensure global cache version exists
        if not hasattr(self, '_cache_version'):
            self._cache_version = 0
        
        # Check if cache exists and is valid
        if hasattr(self, cache_name) and hasattr(self, version_name):
            if getattr(self, version_name) == self._cache_version:
                cached_value = getattr(self, cache_name)
                # For mutable types (set, list, dict), return a copy to prevent cache modification
                if isinstance(cached_value, (set, list, dict)):
                    return type(cached_value)(cached_value)  # Return a copy
                return cached_value
        
        # Compute and cache
        result = func(self)
        setattr(self, cache_name, result)
        setattr(self, version_name, self._cache_version)
        
        # For mutable types (set, list, dict), return a copy to prevent cache modification
        if isinstance(result, (set, list, dict)):
            return type(result)(result)  # Return a copy
        return result
    
    return wrapper


def clears_cache(func: Callable[..., T]) -> Callable[..., T]:
    """
    Decorator to clear all cached properties when a mutation method is called.
    
    This decorator should be applied to methods that modify the object's state
    in a way that would affect cached property values. It increments the global
    cache version, which causes all cached properties to be recomputed on
    their next access.
    
    Parameters
    ----------
    func : Callable[..., T]
        The mutation method that should clear caches.
    
    Returns
    -------
    Callable[..., T]
        Wrapped function that clears caches after execution.
    
    Examples
    --------
    >>> class MyClass:
    ...     def __init__(self):
    ...         self._data = []
    ...         self._cache_version = 0
    ...
    ...     @cached_property
    ...     def count(self):
    ...         return len(self._data)
    ...
    ...     @clears_cache
    ...     def add_item(self, item):
    ...         self._data.append(item)
    """
    @wraps(func)
    def wrapper(self: Any, *args: Any, **kwargs: Any) -> T:
        """
        Wrapper that clears caches after method execution.
        
        Parameters
        ----------
        self : Any
            The instance.
        *args : Any
            Positional arguments.
        **kwargs : Any
            Keyword arguments.
        
        Returns
        -------
        T
            The result of the original method.
        """
        # Execute the mutation method
        result = func(self, *args, **kwargs)
        
        # Increment global version to clear all caches
        if hasattr(self, '_cache_version'):
            self._cache_version += 1
        else:
            self._cache_version = 1
        
        return result
    
    return wrapper


def clear_all_caches(obj: Any) -> None:
    """
    Clear all cached properties on an object.
    
    This should be called from validation methods to ensure caches
    are cleared when state changes are detected.
    
    Parameters
    ----------
    obj : Any
        The object whose caches should be cleared.
    
    Examples
    --------
    >>> class MyClass:
    ...     def __init__(self):
    ...         self._cache_version = 0
    ...
    ...     @cached_property
    ...     def expensive(self):
    ...         return sum(self.data)
    ...
    ...     def _validate(self):
    ...         from phylozoo.utils.cache import clear_all_caches
    ...         clear_all_caches(self)  # Clear caches when validating
    """
    if not hasattr(obj, '_cache_version'):
        obj._cache_version = 0
    else:
        obj._cache_version += 1

