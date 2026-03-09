"""
Styling for DirectedMultiGraph plots.

This module provides styling configuration for DirectedMultiGraph visualizations.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class BaseStyle:
    """
    Base styling configuration for all graph and network plots.

    This class provides the common styling options shared across all
    visualization types.

    Examples
    --------
    >>> style = BaseStyle(node_color='blue')
    >>> style.node_color
    'blue'
    
    Attributes
    ----------
    node_color : str
        Color for nodes. Default is 'lightblue'.
    node_size : float
        Size of nodes. Default is 500.0.
    edge_color : str
        Color for edges. Default is 'gray'.
    edge_width : float
        Width of edges. Default is 2.0.
    with_labels : bool
        Whether to show node labels. Default is True.
    label_offset : float
        Offset for labels from nodes. Default is 0.12.
    label_font_size : float
        Font size for labels. Default is 10.0.
    label_color : str
        Color for labels. Default is 'black'.
    """

    node_color: str = 'lightblue'
    node_size: float = 500.0
    node_edge_color: str = 'black'
    node_edge_width: float = 1.5
    edge_color: str = 'gray'
    edge_width: float = 2.0
    with_labels: bool = True
    label_offset: float = 0.12
    label_font_size: float = 10.0
    label_color: str = 'black'


@dataclass
class DMGraphStyle(BaseStyle):
    """
    Styling configuration for DirectedMultiGraph plots.

    This class extends BaseStyle with DirectedMultiGraph-specific options.

    Examples
    --------
    >>> style = DMGraphStyle(node_color='blue')
    >>> style.node_color
    'blue'
    
    Attributes
    ----------
    node_color : str
        Color for nodes. Default is 'lightblue'.
    node_size : float
        Size of nodes. Default is 500.0.
    edge_color : str
        Color for edges. Default is 'gray'.
    edge_width : float
        Width of edges. Default is 2.0.
    with_labels : bool
        Whether to show node labels. Default is True.
    label_offset : float
        Offset for labels from nodes. Default is 0.12.
    label_font_size : float
        Font size for labels. Default is 10.0.
    label_color : str
        Color for labels. Default is 'black'.
    """

    pass


def default_style() -> DMGraphStyle:
    """
    Get the default style configuration for DirectedMultiGraph.

    Returns
    -------
    DMGraphStyle
        Default style configuration.

    Examples
    --------
    >>> style = default_style()
    >>> style.node_color
    'lightblue'
    """
    return DMGraphStyle()
