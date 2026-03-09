"""
Styling for DirectedPhyNetwork plots.

This module provides styling configuration for DirectedPhyNetwork visualizations.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..d_multigraph.style import DMGraphStyle


@dataclass
class DNetStyle(DMGraphStyle):
    """
    Styling configuration for DirectedPhyNetwork plots.

    This class extends DMGraphStyle with DirectedPhyNetwork-specific options,
    including support for leaves, hybrid nodes, and hybrid edges.

    Examples
    --------
    >>> style = DNetStyle(node_color='blue', leaf_color='green')
    >>> style.node_color
    'blue'
    """

    leaf_color: str = 'lightblue'
    hybrid_color: str = 'lightblue'
    leaf_size: float | None = None
    hybrid_edge_color: str = 'red'

    def copy(self) -> 'DNetStyle':
        """
        Create a copy of this style.

        Returns
        -------
        DNetStyle
            A new DNetStyle instance with the same values.

        Examples
        --------
        >>> style = DNetStyle(node_color='blue')
        >>> style2 = style.copy()
        >>> style2.node_color
        'blue'
        """
        return DNetStyle(
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


def default_style() -> DNetStyle:
    """
    Get the default style configuration for DirectedPhyNetwork.

    Returns
    -------
    DNetStyle
        Default style configuration.

    Examples
    --------
    >>> style = default_style()
    >>> style.node_color
    'lightblue'
    """
    return DNetStyle()
