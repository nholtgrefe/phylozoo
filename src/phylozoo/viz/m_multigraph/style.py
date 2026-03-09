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
