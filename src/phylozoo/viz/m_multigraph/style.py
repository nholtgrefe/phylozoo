"""
Styling for MixedMultiGraph plots.

This module provides styling configuration for MixedMultiGraph visualizations.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..d_multigraph.style import BaseStyle


@dataclass
class MGraphStyle(BaseStyle):
    """
    Styling configuration for MixedMultiGraph plots.

    This class extends BaseStyle with MixedMultiGraph-specific options.
    MixedMultiGraphs can have both directed and undirected edges.

    Examples
    --------
    >>> style = MGraphStyle(node_color='blue')
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
        Offset for labels from nodes. Default is 0.1.
    label_font_size : float
        Font size for labels. Default is 10.0.
    label_color : str
        Color for labels. Default is 'black'.
    """

    pass


def default_style() -> MGraphStyle:
    """
    Get the default style configuration for MixedMultiGraph.

    Returns
    -------
    MGraphStyle
        Default style configuration.

    Examples
    --------
    >>> style = default_style()
    >>> style.node_color
    'lightblue'
    """
    return MGraphStyle()
