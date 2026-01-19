"""
Public API for SemiDirectedPhyNetwork plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import Circle

from phylozoo.utils.exceptions import PhyloZooLayoutError

from .style import SDNetStyle, default_style
from phylozoo.viz._types import EdgeRoute
from .layout import compute_nx_layout, compute_pz_radial_layout

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork


def plot_sdnetwork(
    network: 'SemiDirectedPhyNetwork',
    layout: str = 'twopi',
    style: SDNetStyle | None = None,
    ax: Any | None = None,
    show: bool = False,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a SemiDirectedPhyNetwork.

    This is the main public API function for plotting semi-directed networks.
    It handles layout computation, styling, and rendering using matplotlib.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The network to plot.
    layout : str, optional
        Layout algorithm. Custom PhyloZoo layouts: 'pz-radial', 'pz-planar'.
        NetworkX layouts: 'spring', 'circular', 'kamada_kawai', etc.
        Graphviz layouts: 'dot', 'neato', 'twopi', etc.
        By default 'twopi'.
    style : NetworkStyle, optional
        Styling configuration. If None, uses default style.
        By default None.
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on. If None, creates new figure and axes.
        By default None.
    show : bool, optional
        If True, automatically display the plot using plt.show().
        By default False.
    **layout_kwargs
        Additional parameters for layout computation.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz import plot_sdnetwork
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> ax = plot_sdnetwork(net)
    """
    # Get or create style
    if style is None:
        style = default_style()

    # Create or use existing axes
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Compute layout
    if layout == 'pz-radial':
        computed_layout = compute_pz_radial_layout(network, **layout_kwargs)
    elif layout.startswith('pz-'):
        raise PhyloZooLayoutError(
            f"Unknown PhyloZoo layout: '{layout}'. "
            "Supported PhyloZoo layouts: 'pz-radial'"
        )
    else:
        # NetworkX or Graphviz layout
        computed_layout = compute_nx_layout(network, layout=layout, **layout_kwargs)

    # Render edges
    network_obj = computed_layout.network
    positions = computed_layout.positions

    # Render all edges
    for (u, v, key), route in computed_layout.edge_routes.items():
        _draw_edge(ax, route, style, computed_layout)

    # Determine node types

    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes

    # Render nodes
    for node, position in positions.items():
        if node in leaves:
            node_type = 'leaf'
        elif node in hybrid_nodes:
            node_type = 'hybrid'
        else:
            node_type = 'tree'

        _draw_node(ax, node, position, node_type, style)

        # Add labels if enabled
        if style.with_labels:
            label = network_obj.get_label(node)
            if label:
                # For radial layout, position labels radially outward for leaves
                if layout == 'pz-radial' and node_type == 'leaf':
                    _draw_label_radial(ax, position, label, style)
                else:
                    _draw_label(ax, position, label, style, node_type)

    # Set axis properties
    ax.set_aspect('equal')
    ax.axis('off')

    # Show plot if requested
    if show:
        plt.show()

    return ax


def _draw_edge(
    ax: Any,
    route: EdgeRoute,
    style: SDNetStyle,
    layout: Any,
) -> None:
    """Draw a single edge."""
    points = route.points
    if not points:
        return

    # Determine edge color
    edge_color = (
        style.hybrid_edge_color
        if route.edge_type.is_hybrid
        else style.edge_color
    )

    # Draw edge
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.plot(
        xs,
        ys,
        color=edge_color,
        linewidth=style.edge_width,
        zorder=1,
    )

    # Add arrow for directed edges only
    if route.edge_type.is_directed and len(points) >= 2:
        ax.annotate(
            '',
            xy=(xs[-1], ys[-1]),
            xytext=(xs[-2], ys[-2]),
            arrowprops=dict(
                arrowstyle='->',
                color=edge_color,
                lw=style.edge_width,
            ),
            zorder=2,
        )


def _draw_node(
    ax: Any,
    node: Any,
    position: tuple[float, float],
    node_type: str,
    style: SDNetStyle,
) -> Circle:
    """Draw a single node."""
    x, y = position

    # Determine node color and size
    if node_type == 'leaf':
        color = style.leaf_color
        size = style.leaf_size
    elif node_type == 'hybrid':
        color = style.hybrid_color
        size = style.node_size
    else:
        color = style.node_color
        size = style.node_size

    # Convert size to radius
    radius = (size / 1000.0) ** 0.5 * 0.03

    circle = Circle(
        (x, y),
        radius=radius,
        color=color,
        zorder=3,
    )
    ax.add_patch(circle)
    return circle


def _draw_label(
    ax: Any,
    position: tuple[float, float],
    text: str,
    style: SDNetStyle,
    node_type: str = 'leaf',
) -> Any:
    """Add a text label."""
    x, y = position

    # Position labels based on node type
    if node_type == 'leaf':
        # Leaves: place below (negative y offset)
        offset_x = 0.0
        offset_y = -style.label_offset
        ha = 'center'
        va = 'top'
    else:
        # Internal nodes: place to the right
        offset_x = style.label_offset
        offset_y = 0.0
        ha = 'left'
        va = 'center'

    return ax.text(
        x + offset_x,
        y + offset_y,
        text,
        fontsize=style.label_font_size,
        color=style.label_color,
        ha=ha,
        va=va,
        zorder=4,
    )


def _draw_label_radial(
    ax: Any,
    position: tuple[float, float],
    text: str,
    style: SDNetStyle,
) -> Any:
    """Add a text label positioned radially outward (for radial layout)."""
    import math
    x, y = position
    
    # Calculate angle and radius
    angle = math.atan2(y, x)
    radius = math.sqrt(x * x + y * y)
    
    # Position label radially outward
    label_radius = radius + style.label_offset
    label_x = label_radius * math.cos(angle)
    label_y = label_radius * math.sin(angle)
    
    # Determine horizontal alignment based on angle
    if abs(angle) < math.pi / 6 or abs(angle) > 5 * math.pi / 6:
        ha = 'center'
    elif angle > 0:
        ha = 'left'
    else:
        ha = 'right'
    
    # Vertical alignment
    if abs(angle - math.pi / 2) < math.pi / 6:
        va = 'bottom'
    elif abs(angle + math.pi / 2) < math.pi / 6:
        va = 'top'
    else:
        va = 'center'

    return ax.text(
        label_x,
        label_y,
        text,
        fontsize=style.label_font_size,
        color=style.label_color,
        ha=ha,
        va=va,
        zorder=4,
    )
