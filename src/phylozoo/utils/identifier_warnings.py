"""
Warning utilities for Python keyword identifiers and values.

This module provides functions to warn users when they use Python keywords
as identifiers or when None is used as a value.
"""

import keyword
import warnings
from typing import Any

from .exceptions import PhyloZooIdentifierWarning


def warn_on_keyword(value: Any, context: str) -> None:
    """
    Warn if value is a Python keyword (e.g., None, True, False).

    Parameters
    ----------
    value : Any
        Value to check if it's a Python keyword.
    context : str
        Context description (e.g., "Identifier", "Key", "Name").

    Examples
    --------
    >>> warn_on_keyword("for", "Identifier")
    UserWarning: Identifier 'for' is a Python keyword...
    """
    try:
        as_str = str(value)
    except Exception:
        return
    if keyword.iskeyword(as_str):
        warnings.warn(
            f"{context} '{value}' is a Python keyword; using it as an identifier may cause unexpected behavior.",
            PhyloZooIdentifierWarning,
            stacklevel=3,
        )


def warn_on_none_value(value: Any, context: str) -> None:
    """
    Warn if value is None (Python keyword).

    Parameters
    ----------
    value : Any
        Value to check.
    context : str
        Context description (e.g., "Attribute 'weight'", "Parameter 'x'").

    Examples
    --------
    >>> warn_on_none_value(None, "Attribute 'weight'")
    UserWarning: Attribute 'weight' has value None (Python keyword)...
    """
    if value is None:
        warnings.warn(
            f"{context} has value None (Python keyword); using None as a value may cause unexpected behavior.",
            PhyloZooIdentifierWarning,
            stacklevel=3,
        )

