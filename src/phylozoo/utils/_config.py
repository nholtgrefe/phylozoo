"""
Configuration module.

This module provides configuration settings and validation functions.
"""

from typing import Callable, Optional

# Global validation flag
_validate_enabled: bool = True


def set_validate(enabled: bool) -> None:
    """
    Enable or disable validation.

    Parameters
    ----------
    enabled : bool
        If True, enable validation; if False, disable it
    """
    global _validate_enabled
    _validate_enabled = enabled


def validate() -> bool:
    """
    Check if validation is enabled.

    Returns
    -------
    bool
        True if validation is enabled, False otherwise
    """
    return _validate_enabled
