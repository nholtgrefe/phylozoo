"""
Styling for SemiDirectedPhyNetwork plots.

This module provides styling configuration for SemiDirectedPhyNetwork visualizations.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..m_multigraph.style import MGraphStyle


@dataclass
class SDNetStyle(MGraphStyle):
    """
    Styling configuration for SemiDirectedPhyNetwork plots.

    This class extends MGraphStyle with SemiDirectedPhyNetwork-specific options,
    including support for leaves, hybrid nodes, and hybrid edges.

    Examples
    --------
    >>> style = SDNetStyle(node_color='blue', leaf_color='green')
    >>> style.node_color
    'blue'
    
    Attributes
    ----------
    node_color : str
        Color for internal tree nodes. Default is 'lightblue'.
    leaf_color : str
        Color for leaf nodes. Default is 'lightblue' (same as node_color).
    hybrid_color : str
        Color for hybrid nodes. Default is 'lightblue' (same as node_color).
    node_size : float
        Size of internal nodes. Default is 500.0.
    leaf_size : float | None
        Size of leaf nodes. If None, uses node_size. Default is None.
    edge_color : str
        Color for tree edges. Default is 'gray'.
    hybrid_edge_color : str
        Color for hybrid edges. Default is 'red'.
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

    leaf_color: str = 'lightblue'
    hybrid_color: str = 'lightblue'
    leaf_size: float | None = None
    hybrid_edge_color: str = 'red'

    def copy(self) -> 'SDNetStyle':
        """
        Create a copy of this style.

        Returns
        -------
        SDNetStyle
            A new SDNetStyle instance with the same values.

        Examples
        --------
        >>> style = SDNetStyle(node_color='blue')
        >>> style2 = style.copy()
        >>> style2.node_color
        'blue'
        """
        return SDNetStyle(
            node_color=self.node_color,
            leaf_color=self.leaf_color,
            hybrid_color=self.hybrid_color,
            node_size=self.node_size,
            leaf_size=self.leaf_size,
            node_edge_color=self.node_edge_color,
            node_edge_width=self.node_edge_width,
            edge_color=self.edge_color,
            hybrid_edge_color=self.hybrid_edge_color,
            edge_width=self.edge_width,
            with_labels=self.with_labels,
            label_offset=self.label_offset,
            label_font_size=self.label_font_size,
            label_color=self.label_color,
        )


def default_style() -> SDNetStyle:
    """
    Get the default style configuration for SemiDirectedPhyNetwork.

    Returns
    -------
    SDNetStyle
        Default style configuration.

    Examples
    --------
    >>> style = default_style()
    >>> style.node_color
    'lightblue'
    """
    return SDNetStyle()
