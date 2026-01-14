"""
Validation suppression utilities.

This module provides context managers and class decorators to temporarily disable
validation methods in a controlled, nestable, and context-local way.

Suppression can be specified at two levels:
- **Class level**: Suppress validation for specific classes using fnmatch patterns
- **Method level**: Suppress validation for specific methods using fnmatch patterns

Examples
--------
>>> @validation_aware(allowed=["validate", "_check_inner"], default=["validate"])
... class MyClass:
...     def __init__(self) -> None:
...         self.validate_calls = 0
...         self.check_calls = 0
...         # You can also suppress this initial validate via no_validation().
...         self.validate()  # triggers validate (and _check_inner) once on init
...
...     def validate(self) -> bool:
...         self.validate_calls += 1
...         self._check_inner()
...         return True
...
...     def _check_inner(self) -> bool:
...         self.check_calls += 1
...         return True
...
>>> obj = MyClass()
>>> (obj.validate_calls, obj.check_calls)
(1, 1)
>>> # Suppress validation during construction (including the init call)
>>> # Uses class defaults when no methods specified
>>> with no_validation():
...     obj2 = MyClass()
>>> (obj2.validate_calls, obj2.check_calls)
(0, 0)
>>> # Suppress only the inner check; outer validate still runs
>>> with no_validation(methods=["_check_inner"]):
...     obj.validate()
>>> (obj.validate_calls, obj.check_calls)
(2, 1)
>>> # Suppress both validate and its inner call
>>> with no_validation(methods=["validate"]):
...     obj.validate()
>>> (obj.validate_calls, obj.check_calls)
(2, 1)
>>> # Suppress validation for specific classes only
>>> @validation_aware(allowed=["validate"], default=["validate"])
... class ClassA:
...     def validate(self) -> bool:
...         return True
>>> @validation_aware(allowed=["validate"], default=["validate"])
... class ClassB:
...     def validate(self) -> bool:
...         return True
>>> with no_validation(classes=["ClassA"]):
...     a = ClassA()  # validate suppressed
...     b = ClassB()  # validate NOT suppressed
>>> # Suppress specific methods for specific classes
>>> with no_validation(classes=["ClassA"], methods=["validate"]):
...     a = ClassA()  # validate suppressed
>>> # Use wildcards for classes
>>> with no_validation(classes=["*Generator"]):  # All classes ending in "Generator"
...     gen = DirectedGenerator(...)  # validate suppressed
>>> # Nested suppression: outer block suppresses validate; inner would suppress check,
>>> # but validate itself is already skipped
>>> with no_validation(methods=["validate"]):
...     with no_validation(methods=["_check_inner"]):
...         obj.validate()
>>> (obj.validate_calls, obj.check_calls)
(2, 1)
>>> # Outside the context, validation runs normally
>>> obj.validate()
True
>>> (obj.validate_calls, obj.check_calls)
(3, 2)
"""

from __future__ import annotations

import contextlib
import contextvars
import fnmatch
from functools import wraps
from typing import Callable, Iterable, Iterator, List, Sequence, Type, TypeVar

from phylozoo.utils.exceptions import PhyloZooValueError

T = TypeVar("T")

_validation_disabled: contextvars.ContextVar[bool] = contextvars.ContextVar(
    "_validation_disabled", default=False
)
_suppression_stack: contextvars.ContextVar[list[tuple[list[str] | None, list[str] | None]]] = contextvars.ContextVar(
    "_suppression_stack", default=[]
)


def _active_patterns(default_patterns: Sequence[str]) -> tuple[list[str] | None, list[str] | None]:
    """
    Return the active suppression patterns for the current context.

    Parameters
    ----------
    default_patterns : Sequence[str]
        Default patterns defined on the class.

    Returns
    -------
    tuple[list[str] | None, list[str] | None]
        Tuple of (class_patterns, method_patterns) currently active for suppression.
        If None, means all classes/methods (or use defaults for methods).
    """
    stack = _suppression_stack.get()
    if stack:
        return stack[-1]
    return (None, list(default_patterns))


@contextlib.contextmanager
def no_validation(
    classes: Iterable[str] | None = None,
    methods: Iterable[str] | None = None,
) -> Iterator[None]:
    """
    Temporarily disable validation within the current context.

    Parameters
    ----------
    classes : Iterable[str], optional
        Class names or fnmatch patterns to suppress validation for.
        If None, considers any validation_aware class.
        Examples: ``["DirectedGenerator"]``, ``["*Generator"]``, ``["Directed*"]``
    methods : Iterable[str], optional
        Method names or fnmatch patterns to suppress.
        If None, uses the class defaults supplied by ``validation_aware``
        (active defaults in this context).
        Examples: ``["validate"]``, ``["_validate_*"]``, ``["validate", "_check"]``

    Yields
    ------
    None
        Context scope where validation is suppressed for matching classes and methods.

    Examples
    --------
    >>> # Suppress all validation (uses class defaults)
    >>> with no_validation():
    ...     obj = MyClass()
    
    >>> # Suppress only specific methods for all classes
    >>> with no_validation(methods=["validate"]):
    ...     obj.validate()
    
    >>> # Suppress only specific classes (uses their defaults)
    >>> with no_validation(classes=["DirectedGenerator"]):
    ...     gen = DirectedGenerator(...)  # validate suppressed
    ...     net = DirectedPhyNetwork(...)  # validate NOT suppressed
    
    >>> # Suppress specific methods for specific classes
    >>> with no_validation(classes=["DirectedGenerator"], methods=["validate"]):
    ...     gen = DirectedGenerator(...)  # validate suppressed
    """
    previous_disabled = _validation_disabled.get()
    stack_snapshot = list(_suppression_stack.get())
    
    # Get current patterns from stack if available
    if stack_snapshot:
        current_classes, current_methods = stack_snapshot[-1]
    else:
        current_classes, current_methods = None, None
    
    # Determine new patterns
    new_classes = list(classes) if classes is not None else current_classes
    new_methods = list(methods) if methods is not None else current_methods

    stack_snapshot.append((new_classes, new_methods))
    _suppression_stack.set(stack_snapshot)
    _validation_disabled.set(True)
    try:
        yield
    finally:
        stack_snapshot = _suppression_stack.get()
        if stack_snapshot:
            stack_snapshot.pop()
        _suppression_stack.set(stack_snapshot)
        _validation_disabled.set(previous_disabled)


def validation_aware(
    allowed: Sequence[str],
    default: Sequence[str] | None = None,
) -> Callable[[Type[T]], Type[T]]:
    """
    Decorate a class so selected methods honor validation suppression.

    Parameters
    ----------
    allowed : Sequence[str]
        Method names or fnmatch patterns that may be suppressed.
    default : Sequence[str], optional
        Patterns suppressed by default for this class (must match at least one
        allowed pattern).

    Returns
    -------
    Callable[[Type[T]], Type[T]]
        Class decorator that wraps allowed methods.

    Raises
    ------
    PhyloZooValueError
        If a default pattern does not match any allowed pattern.
    """
    allowed_list = list(allowed)
    default_list = list(default or [])

    invalid_defaults = [
        pattern
        for pattern in default_list
        if not any(fnmatch.fnmatch(pattern, allowed_pattern) for allowed_pattern in allowed_list)
    ]
    if invalid_defaults:
        raise PhyloZooValueError(
            "Default suppression patterns must be a subset of allowed patterns; "
            f"invalid: {invalid_defaults}"
        )

    def decorator(cls: Type[T]) -> Type[T]:
        cls._validation_allowed = allowed_list  # type: ignore[attr-defined]
        cls._validation_default = default_list  # type: ignore[attr-defined]

        for name in dir(cls):
            if not any(fnmatch.fnmatch(name, pattern) for pattern in allowed_list):
                continue

            method = getattr(cls, name)
            if not callable(method):
                continue

            @wraps(method)
            def wrapper(self: T, *args: object, __m: Callable = method, __name: str = name, **kwargs: object):
                if _validation_disabled.get():
                    # Get class name for class pattern matching
                    class_name = self.__class__.__name__
                    
                    # Get active patterns (class_patterns, method_patterns)
                    class_patterns, method_patterns = _active_patterns(
                        getattr(self, "_validation_default", default_list)
                    )
                    
                    # Check if class matches (if class patterns specified)
                    if class_patterns is not None:
                        class_matches = any(
                            fnmatch.fnmatch(class_name, pattern) for pattern in class_patterns
                        )
                        if not class_matches:
                            # Class doesn't match, don't suppress
                            return __m(self, *args, **kwargs)
                    
                    # Check if method matches
                    if method_patterns is not None:
                        # Use specified method patterns
                        method_matches = any(
                            fnmatch.fnmatch(__name, pattern) for pattern in method_patterns
                        )
                        if not method_matches:
                            # Method doesn't match, don't suppress
                            return __m(self, *args, **kwargs)
                    else:
                        # Use class defaults - check if method matches any default pattern
                        class_defaults = getattr(self, "_validation_default", default_list)
                        if class_defaults:
                            method_matches = any(
                                fnmatch.fnmatch(__name, pattern) for pattern in class_defaults
                            )
                            if not method_matches:
                                # Method doesn't match defaults, don't suppress
                                return __m(self, *args, **kwargs)
                        else:
                            # No method patterns specified and no class defaults - don't suppress
                            return __m(self, *args, **kwargs)
                    
                    # Check if method is in allowed list
                    if not any(fnmatch.fnmatch(__name, allowed_pattern) for allowed_pattern in allowed_list):
                        return __m(self, *args, **kwargs)
                    
                    # All checks passed - suppress validation
                    return None
                return __m(self, *args, **kwargs)

            setattr(cls, name, wrapper)

        return cls

    return decorator

