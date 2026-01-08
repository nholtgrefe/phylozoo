"""
Default style configurations for DirectedPhyNetwork plots.

This module provides default style configurations.
"""

from .style import NetworkStyle


def default_style() -> NetworkStyle:
    """
    Get the default style configuration.

    Returns
    -------
    NetworkStyle
        Default style configuration.

    Examples
    --------
    >>> style = default_style()
    >>> style.node_color
    'lightblue'
    """
    return NetworkStyle()
