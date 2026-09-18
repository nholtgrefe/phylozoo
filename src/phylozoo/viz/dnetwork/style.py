"""
Styling for DirectedPhyNetwork plots.

This module provides styling configuration for DirectedPhyNetwork visualizations.
The defaults match those of semi-directed network plots (white tree nodes, black
leaves, pink hybrid nodes, red hybrid edges, labels aligned with their edge).
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

    Attributes
    ----------
    arrows : str | None
        Which edges get an arrowhead: ``'all'``, ``'hybrid'`` (only edges into
        hybrid nodes) or ``'none'``. ``None`` (default) chooses ``'hybrid'`` for
        the layered ``pz-cladogram`` layout, where the direction is implied by the
        drawing, and ``'all'`` for other layouts.
    label_rotation : float | None
        ``None`` (default) aligns each label with the direction of its edge;
        a float applies that fixed rotation (in degrees) to all labels.

    Examples
    --------
    >>> style = DNetStyle(node_color='blue', leaf_color='green')
    >>> style.node_color
    'blue'
    """

    node_color: str = "white"
    leaf_color: str = "#0a0a0a"
    hybrid_color: str = "#fcc0bc"
    leaf_size: float | None = 100.0
    hybrid_edge_color: str = "red"
    node_size: float = 80.0
    arrow_head_size: float = 12.0
    label_offset: float = 0.015
    label_rotation: float | None = None
    arrows: str | None = None

    def copy(self) -> "DNetStyle":
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
            arrow_head_size=self.arrow_head_size,
            with_labels=self.with_labels,
            label_offset=self.label_offset,
            label_font_size=self.label_font_size,
            label_color=self.label_color,
            label_rotation=self.label_rotation,
            arrows=self.arrows,
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
    'white'
    """
    return DNetStyle()
