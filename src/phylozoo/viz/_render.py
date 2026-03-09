"""
Shared rendering utilities for PhyloZoo visualization.

This module provides backend-agnostic drawing primitives (edges, nodes, labels)
used by all plot functions. Styles must have BaseStyle attributes; network
styles may add leaf_color, hybrid_color, leaf_size, hybrid_edge_color.

Layout algorithms only determine node positions; node sizes come from the style.
"""

from __future__ import annotations

from typing import Any, Callable, Protocol

import matplotlib.patches as mpatches
from matplotlib.path import Path as MPath
from matplotlib.patches import Circle

from ._types import EdgeRoute


class RenderStyle(Protocol):
    """Protocol for style objects used in rendering."""

    node_color: str
    node_size: float
    edge_color: str
    edge_width: float
    with_labels: bool
    label_offset: float
    label_font_size: float
    label_color: str


# Parallel edge offset constants (shared across all plotters)
BASE_PARALLEL_OFFSET = 0.15
PARALLEL_OFFSET_STEP = 0.1


def _get_edge_color(route: EdgeRoute, style: RenderStyle) -> str:
    """Get edge color from route type and style."""
    if route.edge_type.is_hybrid and hasattr(style, 'hybrid_edge_color'):
        return style.hybrid_edge_color
    return style.edge_color


def _compute_parallel_offset(
    route: EdgeRoute,
    parallel_groups: dict[tuple[Any, Any], list[int]],
    u: Any,
    v: Any,
    key: int,
) -> float:
    """Compute offset for parallel edges to avoid overlap."""
    if not route.edge_type.is_parallel:
        return 0.0
    pair = (u, v) if route.edge_type.is_directed else (min(u, v), max(u, v))
    keys = parallel_groups.get(pair, [])
    if len(keys) <= 1 or key not in keys:
        return 0.0
    parallel_keys = sorted(keys)
    key_index = parallel_keys.index(key)
    direction = 1 if key_index % 2 == 0 else -1
    magnitude = BASE_PARALLEL_OFFSET + (key_index // 2) * PARALLEL_OFFSET_STEP
    return direction * magnitude


def draw_edge(
    ax: Any,
    route: EdgeRoute,
    style: RenderStyle,
    parallel_groups: dict[tuple[Any, Any], list[int]],
    edge_key: tuple[Any, Any, int],
) -> None:
    """
    Draw a single edge.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    route : EdgeRoute
        Edge routing information.
    style : RenderStyle
        Styling configuration. Must have edge_color, edge_width. Network
        styles may have hybrid_edge_color for hybrid edges.
    parallel_groups : dict
        Mapping of (u, v) to list of edge keys for parallel edge offset.
    edge_key : tuple
        (u, v, key) for this edge.
    """
    points = route.points
    if not points:
        return

    edge_color = _get_edge_color(route, style)
    u, v, key = edge_key
    parallel_offset = _compute_parallel_offset(route, parallel_groups, u, v, key)

    if route.edge_type.is_parallel and parallel_offset != 0.0:
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
        if route.edge_type.is_directed:
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
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        ax.plot(
            xs,
            ys,
            color=edge_color,
            linewidth=style.edge_width,
            zorder=1,
        )
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


def draw_node(
    ax: Any,
    position: tuple[float, float],
    node_type: str,
    style: RenderStyle,
) -> Circle:
    """
    Draw a single node.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    position : tuple[float, float]
        (x, y) position.
    node_type : str
        One of 'root', 'leaf', 'tree', 'hybrid', 'generic'. For network
        plots use root/leaf/tree/hybrid; for graph plots use 'generic'.
    style : RenderStyle
        Styling. Network styles may have leaf_color, leaf_size, hybrid_color.
        May have node_edge_color, node_edge_width for border.

    Returns
    -------
    Circle
        The drawn circle patch.
    """
    x, y = position
    color = _get_node_color(node_type, style)
    size = _get_node_size(node_type, style)
    # Linear scaling: size 500 -> ~0.02 radius; size 20000 -> ~0.5 (capped)
    radius = min(size / 1000.0, 12.0) * 0.042

    edgecolor = getattr(style, 'node_edge_color', 'black')
    linewidth = getattr(style, 'node_edge_width', 1.5)

    circle = Circle(
        (x, y),
        radius=radius,
        facecolor=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=3,
    )
    ax.add_patch(circle)
    return circle


def _get_node_color(node_type: str, style: RenderStyle) -> str:
    """Get node color from node type and style."""
    if node_type == 'leaf' and hasattr(style, 'leaf_color'):
        return style.leaf_color
    if node_type == 'hybrid' and hasattr(style, 'hybrid_color'):
        return style.hybrid_color
    return style.node_color


def _get_node_size(node_type: str, style: RenderStyle) -> float:
    """Get node size from node type and style.

    For leaf nodes, uses leaf_size if set; otherwise falls back to node_size.
    """
    if node_type == 'leaf' and hasattr(style, 'leaf_size') and style.leaf_size is not None:
        return style.leaf_size
    return style.node_size


def draw_label(
    ax: Any,
    position: tuple[float, float],
    text: str,
    style: RenderStyle,
    node_type: str = 'leaf',
    center: tuple[float, float] | None = None,
) -> Any:
    """
    Add a text label, placed outward from the layout center for consistent placement.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    position : tuple[float, float]
        (x, y) position.
    text : str
        Label text.
    style : RenderStyle
        Styling configuration.
    node_type : str, optional
        Used for fallback when node coincides with center. By default 'leaf'.
    center : tuple[float, float] | None, optional
        Layout center for outward placement. If None, uses origin (0, 0).
        By default None.

    Returns
    -------
    matplotlib.text.Text
        The text object.
    """
    import math

    x, y = position
    cx, cy = center if center is not None else (0.0, 0.0)
    dx = x - cx
    dy = y - cy
    dist = math.sqrt(dx * dx + dy * dy)

    if dist > 1e-6:
        # Place label outward from center
        scale = style.label_offset / dist
        offset_x = dx * scale
        offset_y = dy * scale
        if abs(dx) < 1e-6:
            ha = 'center'
        elif dx > 0:
            ha = 'left'
        else:
            ha = 'right'
        if abs(dy) < 1e-6:
            va = 'center'
        elif dy > 0:
            va = 'bottom'
        else:
            va = 'top'
    else:
        offset_x = 0.0
        offset_y = -style.label_offset
        ha = 'center'
        va = 'top'

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


def draw_label_radial(
    ax: Any,
    position: tuple[float, float],
    text: str,
    style: RenderStyle,
) -> Any:
    """
    Add a text label positioned radially outward (for radial layout).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    position : tuple[float, float]
        (x, y) position of node.
    text : str
        Label text.
    style : RenderStyle
        Styling configuration.

    Returns
    -------
    matplotlib.text.Text
        The text object.
    """
    import math

    x, y = position
    angle = math.atan2(y, x)
    radius = math.sqrt(x * x + y * y)
    label_radius = radius + style.label_offset
    label_x = label_radius * math.cos(angle)
    label_y = label_radius * math.sin(angle)

    if abs(angle) < math.pi / 6 or abs(angle) > 5 * math.pi / 6:
        ha = 'center'
    elif angle > 0:
        ha = 'left'
    else:
        ha = 'right'

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


def render_layout(
    ax: Any,
    edge_routes: dict[tuple[Any, Any, int], EdgeRoute],
    positions: dict[Any, tuple[float, float]],
    style: RenderStyle,
    center: tuple[float, float],
    get_node_type: Callable[[Any], str],
    get_label: Callable[[Any], str | None],
    radial_labels_for_leaves: bool = False,
) -> None:
    """
    Shared render loop: draw edges, nodes, and labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    edge_routes : dict
        (u, v, key) -> EdgeRoute.
    positions : dict
        Node ID -> (x, y).
    style : RenderStyle
        Styling configuration.
    center : tuple[float, float]
        Layout center for label placement.
    get_node_type : callable
        node_id -> 'root'|'leaf'|'tree'|'hybrid'|'generic'.
    get_label : callable
        node_id -> label string or None.
    radial_labels_for_leaves : bool, optional
        If True, use draw_label_radial for leaf nodes. By default False.
    """
    parallel_groups = build_parallel_groups(edge_routes)

    for (u, v, key), route in edge_routes.items():
        draw_edge(ax, route, style, parallel_groups, (u, v, key))

    for node, position in positions.items():
        node_type = get_node_type(node)
        draw_node(ax, position, node_type, style)
        if style.with_labels:
            label = get_label(node)
            if label:
                if radial_labels_for_leaves and node_type == 'leaf':
                    draw_label_radial(ax, position, label, style)
                else:
                    draw_label(ax, position, label, style, node_type, center=center)

    ax.set_aspect('equal')
    ax.axis('off')


def build_parallel_groups(
    edge_routes: dict[tuple[Any, Any, int], EdgeRoute],
) -> dict[tuple[Any, Any], list[int]]:
    """
    Build mapping of edge pair to keys for parallel edge offset calculation.

    For directed edges uses (u, v) as key; for undirected uses (min(u,v), max(u,v)).

    Parameters
    ----------
    edge_routes : dict
        Edge routing mapping (u, v, key) -> EdgeRoute.

    Returns
    -------
    dict
        Mapping of (u,v) or (min,max) to list of edge keys.
    """
    groups: dict[tuple[Any, Any], list[int]] = {}
    for (u, v, key), route in edge_routes.items():
        if route.edge_type.is_parallel:
            pair = (u, v) if route.edge_type.is_directed else (min(u, v), max(u, v))
            if pair not in groups:
                groups[pair] = []
            groups[pair].append(key)
    return groups
