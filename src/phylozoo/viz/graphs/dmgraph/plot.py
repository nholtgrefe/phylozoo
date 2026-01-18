"""
Public API for DirectedMultiGraph plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import Circle

from phylozoo.utils.exceptions import PhyloZooLayoutError

from .layout.nx import compute_nx_layout
from .style import DMGraphStyle, default_style
from phylozoo.viz._types import EdgeRoute

if TYPE_CHECKING:
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


def plot_dmgraph(
    graph: 'DirectedMultiGraph',
    layout: str = 'spring',
    style: DMGraphStyle | None = None,
    ax: Any | None = None,
    show: bool = False,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a DirectedMultiGraph.

    This is the main public API function for plotting directed multigraphs.
    It handles layout computation, styling, and rendering using matplotlib.

    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to plot.
    layout : str, optional
        Layout algorithm. NetworkX layouts: 'spring', 'circular', 'kamada_kawai', etc.
        Graphviz layouts: 'dot', 'neato', 'twopi', etc.
        By default 'spring'.
    style : DMGraphStyle, optional
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
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.viz.graphs.dmgraph import plot_dmgraph
    >>>
    >>> G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
    >>> ax = plot_dmgraph(G, layout='circular')
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
    computed_layout = compute_nx_layout(graph, layout=layout, **layout_kwargs)

    # Render edges
    positions = computed_layout.positions

    # Group parallel edges for offset calculation
    base_parallel_offset = 0.15
    parallel_offset_step = 0.1
    parallel_groups: dict[tuple[Any, Any], list[int]] = {}
    for (u, v, key), route in computed_layout.edge_routes.items():
        if route.edge_type.is_parallel:
            edge_key = (u, v)
            if edge_key not in parallel_groups:
                parallel_groups[edge_key] = []
            parallel_groups[edge_key].append(key)

    # Render all edges
    for (u, v, key), route in computed_layout.edge_routes.items():
        _draw_edge(ax, route, style, parallel_groups, key)

    # Render nodes
    # For graphs, all nodes are treated the same
    for node, position in positions.items():
        _draw_node(ax, node, position, style)

        # Add labels if enabled
        if style.with_labels:
            # Try to get label from graph node attributes
            label = None
            if node in graph._graph.nodes:
                node_attrs = graph._graph.nodes[node]
                label = node_attrs.get('label', str(node))
            
            if label:
                _draw_label(ax, position, label, style)

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
    style: DMGraphStyle,
    parallel_groups: dict[tuple[Any, Any], list[int]],
    key: int,
) -> None:
    """Draw a single edge."""
    points = route.points
    if not points:
        return

    edge_color = style.edge_color

    # Calculate parallel edge offset if needed
    base_parallel_offset = 0.15
    parallel_offset_step = 0.1
    parallel_offset = 0.0
    if route.edge_type.is_parallel:
        # Find the edge key in parallel groups
        for (u_node, v_node), keys in parallel_groups.items():
            if len(keys) > 1 and key in keys:
                parallel_keys = sorted(keys)
                key_index = parallel_keys.index(key)
                direction = 1 if key_index % 2 == 0 else -1
                magnitude = base_parallel_offset + (key_index // 2) * parallel_offset_step
                parallel_offset = direction * magnitude
                break

    # Draw edge
    if route.edge_type.is_parallel and parallel_offset != 0.0:
        # Bezier curve for parallel edges
        px, py = points[0]
        ex, ey = points[-1]
        mid_x = (px + ex) / 2
        mid_y = (py + ey) / 2

        if abs(ex - px) < abs(ey - py):
            cx = mid_x + parallel_offset
            cy = mid_y
        else:
            cx = mid_x
            cy = mid_y + parallel_offset

        verts = [(px, py), (cx, cy), (ex, ey)]
        codes = [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]
        path = MPath(verts, codes)
        patch = mpatches.PathPatch(
            path,
            edgecolor=edge_color,
            facecolor='none',
            linewidth=style.edge_width,
            zorder=1,
        )
        ax.add_patch(patch)

        # Add arrow
        dx = ex - cx
        dy = ey - cy
        ax.annotate(
            '',
            xy=(ex, ey),
            xytext=(ex - 0.1 * dx, ey - 0.1 * dy),
            arrowprops=dict(
                arrowstyle='->',
                color=edge_color,
                lw=style.edge_width,
            ),
            zorder=2,
        )
    else:
        # Straight line
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        ax.plot(
            xs,
            ys,
            color=edge_color,
            linewidth=style.edge_width,
            zorder=1,
        )

        # Add arrow for directed edges
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
    style: DMGraphStyle,
) -> Circle:
    """Draw a single node."""
    x, y = position

    # Convert size to radius
    radius = (style.node_size / 1000.0) ** 0.5 * 0.03

    circle = Circle(
        (x, y),
        radius=radius,
        color=style.node_color,
        zorder=3,
    )
    ax.add_patch(circle)
    return circle


def _draw_label(
    ax: Any,
    position: tuple[float, float],
    text: str,
    style: DMGraphStyle,
) -> Any:
    """Add a text label."""
    x, y = position

    # Position labels to the right
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
