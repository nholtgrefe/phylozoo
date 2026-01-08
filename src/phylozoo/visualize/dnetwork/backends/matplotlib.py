"""
Matplotlib backend for DirectedPhyNetwork plotting.

This module provides a matplotlib implementation of the Backend interface.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import Circle

from ...utils.types import EdgeRoute
from .base import Backend, register_backend

if TYPE_CHECKING:
    from ..layout.base import DAGLayout
    from ..styling.style import NetworkStyle


class MatplotlibBackend(Backend):
    """
    Matplotlib backend for network plotting.

    This backend renders network layouts using matplotlib, producing
    static figures that can be displayed or saved.

    Examples
    --------
    >>> from phylozoo.visualize.dnetwork.backends import MatplotlibBackend
    >>> backend = MatplotlibBackend()
    >>> fig = backend.create_figure()
    >>> ax = backend.create_axes(fig)
    """

    def create_figure(self, **kwargs: Any) -> Any:
        """
        Create a new matplotlib figure.

        Parameters
        ----------
        **kwargs
            Arguments passed to plt.figure().

        Returns
        -------
        matplotlib.figure.Figure
            The created figure.
        """
        self._figure = plt.figure(**kwargs)
        return self._figure

    def create_axes(self, figure: Any | None = None, **kwargs: Any) -> Any:
        """
        Create axes on a figure.

        Parameters
        ----------
        figure : matplotlib.figure.Figure, optional
            Existing figure to create axes on. If None, creates a new figure.
        **kwargs
            Arguments passed to figure.add_subplot() or plt.subplots().

        Returns
        -------
        matplotlib.axes.Axes
            The created axes.
        """
        if figure is not None:
            self._figure = figure
        elif self._figure is None:
            self.create_figure()

        if hasattr(self._figure, 'add_subplot'):
            self._axes = self._figure.add_subplot(111, **kwargs)
        else:
            # Fallback if figure doesn't support add_subplot
            _, self._axes = plt.subplots(**kwargs)
            self._figure = self._axes.figure

        return self._axes

    def render_edge(
        self,
        route: EdgeRoute,
        style: 'NetworkStyle',
        layout: 'DAGLayout',
        parallel_offset: float = 0.0,
    ) -> dict[str, Any]:
        """
        Render a single edge using matplotlib.

        Parameters
        ----------
        route : EdgeRoute
            Edge routing information.
        style : NetworkStyle
            Styling configuration.
        layout : DAGLayout
            The layout containing node positions.
        parallel_offset : float, optional
            Offset for parallel edges. By default 0.0.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'patch', 'arrow', or 'line' keys containing
            matplotlib elements.
        """
        if self._axes is None:
            raise ValueError("Axes must be created before rendering")

        ax = self._axes
        points = route.points
        elements: dict[str, Any] = {}

        # Determine edge color
        # Color edges red if they are hybrid edges OR if they go to a hybrid node
        network = layout.network
        # Find target node by matching the endpoint position (with small tolerance for floating point)
        target_pos = points[-1] if points else None
        is_hybrid_target = False
        if target_pos is not None and hasattr(network, 'hybrid_nodes'):
            for node, pos in layout.positions.items():
                # Check if positions match (with tolerance for floating point precision)
                if (abs(pos[0] - target_pos[0]) < 1e-6 
                    and abs(pos[1] - target_pos[1]) < 1e-6 
                    and node in network.hybrid_nodes):
                    is_hybrid_target = True
                    break
        
        edge_color = (
            style.hybrid_edge_color
            if (route.edge_type.is_hybrid or is_hybrid_target)
            else style.edge_color
        )

        # Draw all edges as straight lines (hybrid edges are no longer curved)
        if route.edge_type.is_parallel:
            # For parallel edges, use Bezier curve
            px, py = points[0]
            ex, ey = points[-1]

            # Control point for the curve (perpendicular to edge direction)
            mid_x = (px + ex) / 2
            mid_y = (py + ey) / 2

            # Use parallel_offset if provided, otherwise use default
            offset = parallel_offset if parallel_offset != 0.0 else 0.15

            # Curve perpendicular to edge direction
            if abs(ex - px) < abs(ey - py):
                # More vertical than horizontal
                cx = mid_x + offset
                cy = mid_y
            else:
                # More horizontal than vertical
                cx = mid_x
                cy = mid_y + offset

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
                edgecolor=edge_color,
                facecolor='none',
                linewidth=style.edge_width,
                zorder=1,
            )
            ax.add_patch(patch)
            elements['patch'] = patch

            # Add arrow at the end
            dx = ex - cx
            dy = ey - cy
            arrow = ax.annotate(
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
            elements['arrow'] = arrow
        else:
            # Regular straight edge
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            line = ax.plot(
                xs,
                ys,
                color=edge_color,
                linewidth=style.edge_width,
                zorder=1,
            )[0]
            elements['line'] = line

            # Add arrow at the end
            if len(points) >= 2:
                arrow = ax.annotate(
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
                elements['arrow'] = arrow

        return elements

    def render_node(
        self,
        node: Any,
        position: tuple[float, float],
        node_type: str,
        style: 'NetworkStyle',
    ) -> Circle:
        """
        Render a single node using matplotlib.

        Parameters
        ----------
        node : Any
            Node identifier.
        position : tuple[float, float]
            (x, y) position of the node.
        node_type : str
            Type of node: 'root', 'leaf', 'tree', or 'hybrid'.
        style : NetworkStyle
            Styling configuration.

        Returns
        -------
        matplotlib.patches.Circle
            The rendered node circle.
        """
        if self._axes is None:
            raise ValueError("Axes must be created before rendering")

        ax = self._axes
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

        # Convert size to radius (approximate conversion)
        radius = (size / 1000.0) ** 0.5 * 0.1

        circle = Circle(
            (x, y),
            radius=radius,
            color=color,
            zorder=3,
        )
        ax.add_patch(circle)
        return circle

    def add_label(
        self,
        position: tuple[float, float],
        text: str,
        style: 'NetworkStyle',
        node_type: str = 'leaf',
    ) -> Any:
        """
        Add a text label using matplotlib.

        Parameters
        ----------
        position : tuple[float, float]
            (x, y) position for the label.
        text : str
            Label text.
        style : NetworkStyle
            Styling configuration.
        node_type : str, optional
            Type of node: 'root', 'leaf', 'tree', or 'hybrid'.
            Used to determine label positioning. By default 'leaf'.

        Returns
        -------
        matplotlib.text.Text
            The text label.
        """
        if self._axes is None:
            raise ValueError("Axes must be created before rendering")

        ax = self._axes
        x, y = position

        # Position labels based on node type for top-down layout
        # For leaves: place below the node (negative y offset)
        # For internal nodes: place to the right (positive x offset)
        if node_type == 'leaf':
            # Leaves: below the node
            label_x = x
            label_y = y - style.label_offset * 2  # Larger offset below
        elif node_type == 'root':
            # Root: above the node
            label_x = x
            label_y = y + style.label_offset * 2
        else:
            # Internal nodes (tree, hybrid): to the right
            label_x = x + style.label_offset * 2  # Larger offset to the right
            label_y = y

        text_obj = ax.text(
            label_x,
            label_y,
            text,
            fontsize=style.label_font_size,
            color=style.label_color,
            zorder=4,
            ha='left' if node_type != 'leaf' else 'center',  # Left-align for internal, center for leaves
            va='top' if node_type == 'leaf' else 'center',  # Top-align for leaves, center for others
        )
        return text_obj

    def show(self) -> None:
        """
        Display the matplotlib figure.
        """
        if self._figure is None:
            raise ValueError("Figure must be created before showing")
        plt.show()

    def save(self, path: str, **kwargs: Any) -> None:
        """
        Save the matplotlib figure to a file.

        Parameters
        ----------
        path : str
            File path to save to.
        **kwargs
            Arguments passed to figure.savefig().
        """
        if self._figure is None:
            raise ValueError("Figure must be created before saving")
        self._figure.savefig(path, **kwargs)


# Register the matplotlib backend
register_backend('matplotlib', MatplotlibBackend)

