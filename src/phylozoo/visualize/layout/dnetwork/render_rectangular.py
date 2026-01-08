"""
Rendering functions for rectangular layouts of DirectedPhyNetwork.

This module provides functions to render/draw layout objects. The render
functions only draw the network structure (nodes and edges) without labels,
colors, or styling - those are handled by the public API.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from .base import RectangularDNetLayout


def render_rectangular_dnet_layout(
    layout: 'RectangularDNetLayout',
    ax: 'Axes',
) -> None:
    """
    Render a rectangular DirectedPhyNetwork layout.

    This function draws only the network structure (nodes as points and edges
    as lines/curves). It does not add labels, colors, or styling - those are
    handled by the public plotting API.

    Parameters
    ----------
    layout : RectangularDNetLayout
        The layout to render.
    ax : matplotlib.axes.Axes
        Matplotlib axes to draw on.

    Notes
    -----
    This is a low-level rendering function. For full plotting with labels,
    colors, and styling, use the public API functions in network_plot.py.

    Examples
    --------
    >>> from phylozoo.visualize.layout.dnetwork import (
    ...     compute_rectangular_dnet_layout,
    ...     render_rectangular_dnet_layout
    ... )
    >>> import matplotlib.pyplot as plt
    >>> 
    >>> layout = compute_rectangular_dnet_layout(network)
    >>> fig, ax = plt.subplots()
    >>> render_rectangular_dnet_layout(layout, ax)
    >>> # Add labels, colors, etc. via public API
    """
    network = layout.network
    positions = layout.positions

    # Draw edges
    for (u, v, key), route in layout.edge_routes.items():
        points = route.points

        if route.edge_type.is_hybrid and route.curve_control:
            # Draw curved edge using Bezier curve
            px, py = points[0]
            cx, cy = route.curve_control
            ex, ey = points[-1]

            # Use quadratic Bezier
            verts = [
                (px, py),
                (cx, cy),
                (ex, ey),
            ]
            codes = [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]

            path = MPath(verts, codes)
            patch = mpatches.PathPatch(
                path,
                edgecolor='black',
                facecolor='none',
                linewidth=1.0,
                zorder=1,
            )
            ax.add_patch(patch)

            # Add arrow at the end
            dx = ex - cx
            dy = ey - cy
            ax.annotate(
                '',
                xy=(ex, ey),
                xytext=(ex - 0.1 * dx, ey - 0.1 * dy),
                arrowprops=dict(
                    arrowstyle='->',
                    color='black',
                    lw=1.0,
                ),
                zorder=2,
            )
        else:
            # Draw polyline for tree edges
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            ax.plot(
                xs,
                ys,
                color='black',
                linewidth=1.0,
                zorder=1,
            )
            # Add arrow at the end
            if len(points) >= 2:
                ax.annotate(
                    '',
                    xy=(xs[-1], ys[-1]),
                    xytext=(xs[-2], ys[-2]),
                    arrowprops=dict(
                        arrowstyle='->',
                        color='black',
                        lw=1.0,
                    ),
                    zorder=2,
                )

    # Draw nodes as simple points
    for node, (x, y) in positions.items():
        ax.plot(
            x,
            y,
            marker='o',
            markersize=5,
            color='black',
            zorder=3,
        )

