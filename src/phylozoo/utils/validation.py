"""
Validation suppression utilities.

This module provides context managers and class decorators to temporarily disable
validation methods in a controlled, nestable, and context-local way.

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
>>> with no_validation():
...     obj2 = MyClass()
>>> (obj2.validate_calls, obj2.check_calls)
(0, 0)
>>> # Suppress only the inner check; outer validate still runs
>>> with no_validation(only=["_check_inner"]):
...     obj.validate()
>>> (obj.validate_calls, obj.check_calls)
(2, 1)
>>> # Suppress both validate and its inner call
>>> with no_validation(only=["validate"]):
...     obj.validate()
>>> (obj.validate_calls, obj.check_calls)
(2, 1)
>>> # Nested suppression: outer block suppresses validate; inner would suppress check,
>>> # but validate itself is already skipped
>>> with no_validation(only=["validate"]):
...     with no_validation(only=["_check_inner"]):
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
from typing import Callable, Iterable, Iterator, List, Optional, Sequence, Type, TypeVar


T = TypeVar("T")

_validation_disabled: contextvars.ContextVar[bool] = contextvars.ContextVar(
    "_validation_disabled", default=False
)
_suppression_stack: contextvars.ContextVar[List[List[str]]] = contextvars.ContextVar(
    "_suppression_stack", default=[]
)


def _active_patterns(default_patterns: Sequence[str]) -> List[str]:
    """
    Return the active suppression patterns for the current context.

    Parameters
    ----------
    default_patterns : Sequence[str]
        Default patterns defined on the class.

    Returns
    -------
    List[str]
        Patterns currently active for suppression.
    """
    stack = _suppression_stack.get()
    if stack:
        return stack[-1]
    return list(default_patterns)


@contextlib.contextmanager
def no_validation(only: Optional[Iterable[str]] = None) -> Iterator[None]:
    """
    Temporarily disable validation within the current context.

    Parameters
    ----------
    only : Iterable[str], optional
        Method names or fnmatch patterns to suppress. If None, uses the class
        defaults supplied by ``validation_aware`` (active defaults in this
        context).

    Yields
    ------
    None
        Context scope where validation is suppressed for matching methods.
    """
    previous_disabled = _validation_disabled.get()
    stack_snapshot = list(_suppression_stack.get())
    current_patterns = stack_snapshot[-1] if stack_snapshot else []
    new_patterns = list(only) if only is not None else list(current_patterns)

    stack_snapshot.append(new_patterns)
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
    default: Optional[Sequence[str]] = None,
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
    ValueError
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
        raise ValueError(
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
                    active = _active_patterns(getattr(self, "_validation_default", default_list))
                    if not active:
                        active = getattr(self, "_validation_default", default_list)
                    # Filter patterns to those that apply to this method and are allowed
                    filtered = [
                        pattern
                        for pattern in active
                        if any(fnmatch.fnmatch(__name, allowed_pattern) for allowed_pattern in allowed_list)
                    ]
                    if any(fnmatch.fnmatch(__name, pattern) for pattern in filtered):
                        return None
                return __m(self, *args, **kwargs)

            setattr(cls, name, wrapper)

        return cls

    return decorator

