"""
PyQtGraph backend for DirectedPhyNetwork plotting.

This module provides a PyQtGraph implementation of the Backend interface.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

try:
    import pyqtgraph as pg
    from PyQt5.QtWidgets import QApplication  # type: ignore[import-untyped]
    HAS_PYQTGAPH = True
except ImportError:
    try:
        import pyqtgraph as pg
        from PyQt6.QtWidgets import QApplication  # type: ignore[import-untyped]
        HAS_PYQTGAPH = True
    except ImportError:
        HAS_PYQTGAPH = False
        pg = None
        QApplication = None

from ...utils.types import EdgeRoute
from .base import Backend, register_backend

if TYPE_CHECKING:
    from ..layout.base import DAGLayout
    from ..styling.style import NetworkStyle


class PyQtGraphBackend(Backend):
    """
    PyQtGraph backend for network plotting.

    This backend renders network layouts using PyQtGraph, producing
    interactive plots that can be zoomed and panned.

    Examples
    --------
    >>> from phylozoo.visualize.dnetwork.backends import PyQtGraphBackend
    >>> backend = PyQtGraphBackend()
    >>> win = backend.create_figure()
    >>> plot = backend.create_axes(win)
    """

    def __init__(self) -> None:
        """Initialize PyQtGraph backend."""
        if not HAS_PYQTGAPH:
            raise ImportError(
                "pyqtgraph is not installed. Install it with: pip install pyqtgraph"
            )
        super().__init__()
        self._window: Any | None = None
        self._plot: Any | None = None

    def create_figure(self, **kwargs: Any) -> Any:
        """
        Create a new PyQtGraph window.

        Parameters
        ----------
        **kwargs
            Arguments passed to pg.GraphicsLayoutWidget().

        Returns
        -------
        pg.GraphicsLayoutWidget
            The created window.
        """
        title = kwargs.pop('title', 'Phylogenetic Network')
        show = kwargs.pop('show', True)
        self._window = pg.GraphicsLayoutWidget(show=show, title=title, **kwargs)
        return self._window

    def create_axes(self, figure: Any | None = None, **kwargs: Any) -> Any:
        """
        Create plot on a window.

        Parameters
        ----------
        figure : pg.GraphicsLayoutWidget, optional
            Existing window to create plot on. If None, creates a new window.
        **kwargs
            Arguments passed to window.addPlot().

        Returns
        -------
        pg.PlotItem
            The created plot.
        """
        if figure is not None:
            self._window = figure
        elif self._window is None:
            self.create_figure()

        title = kwargs.pop('title', 'Phylogenetic Network')
        self._plot = self._window.addPlot(title=title, **kwargs)
        return self._plot

    def render_edge(
        self,
        route: EdgeRoute,
        style: 'NetworkStyle',
        layout: 'DAGLayout',
        parallel_offset: float = 0.0,
    ) -> dict[str, Any]:
        """
        Render a single edge using PyQtGraph.

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
            Dictionary with 'line' key containing PyQtGraph line item.
        """
        if self._plot is None:
            raise ValueError("Plot must be created before rendering")

        plot = self._plot
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

        # Convert color name to hex if needed
        color_hex = self._color_to_hex(edge_color)

        # Draw edges
        if route.edge_type.is_parallel and len(points) >= 2:
            # For parallel edges, use curve
            px, py = points[0]
            ex, ey = points[-1]

            # Control point for the curve
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

            # Create Bezier curve points
            import numpy as np
            t = np.linspace(0, 1, 50)
            curve_x = (1 - t) ** 2 * px + 2 * (1 - t) * t * cx + t ** 2 * ex
            curve_y = (1 - t) ** 2 * py + 2 * (1 - t) * t * cy + t ** 2 * ey

            line = plot.plot(
                curve_x,
                curve_y,
                pen={'color': color_hex, 'width': style.edge_width},
            )
            elements['line'] = line
        else:
            # Straight line
            if len(points) >= 2:
                x_coords = [p[0] for p in points]
                y_coords = [p[1] for p in points]

                line = plot.plot(
                    x_coords,
                    y_coords,
                    pen={'color': color_hex, 'width': style.edge_width},
                )
                elements['line'] = line

        return elements

    def render_node(
        self,
        node: Any,
        position: tuple[float, float],
        node_type: str,
        style: 'NetworkStyle',
    ) -> dict[str, Any]:
        """
        Render a single node using PyQtGraph.

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
        dict[str, Any]
            Dictionary with 'point' key containing PyQtGraph scatter plot item.
        """
        if self._plot is None:
            raise ValueError("Plot must be created before rendering")

        plot = self._plot
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

        # Scale down size for PyQtGraph (style sizes are for matplotlib)
        # PyQtGraph symbolSize is typically 5-20, so divide by ~30
        pyqtgraph_size = max(5, size / 30.0)

        color_hex = self._color_to_hex(color)

        point = plot.plot(
            [x],
            [y],
            pen=None,
            symbol='o',
            symbolSize=pyqtgraph_size,
            symbolBrush=color_hex,
        )

        return {'point': point}

    def add_label(
        self,
        position: tuple[float, float],
        text: str,
        style: 'NetworkStyle',
        node_type: str = 'leaf',
    ) -> Any:
        """
        Add a text label using PyQtGraph.

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
        pg.TextItem
            The created text item.
        """
        if self._plot is None:
            raise ValueError("Plot must be created before rendering")

        plot = self._plot
        x, y = position

        # Determine label offset based on node type
        if node_type == 'leaf':
            # Leaf labels below
            offset_x, offset_y = 0.0, -0.2
        elif node_type == 'root':
            # Root labels above
            offset_x, offset_y = 0.0, 0.2
        elif node_type == 'hybrid':
            # Hybrid labels to the right
            offset_x, offset_y = 0.2, 0.0
        else:
            # Internal node labels to the right
            offset_x, offset_y = 0.2, 0.0

        text_item = pg.TextItem(text, color=style.label_color)
        text_item.setPos(x + offset_x, y + offset_y)
        plot.addItem(text_item)

        return text_item

    def show(self) -> None:
        """
        Display the plot window.

        This method ensures a QApplication exists and shows the window.
        """
        if self._window is None:
            raise ValueError("Window must be created before showing")

        # Ensure QApplication exists
        app = QApplication.instance()
        if app is None:
            app = QApplication([])

        self._window.show()

    def save(self, filename: str, **kwargs: Any) -> None:
        """
        Save the plot to a file.

        Parameters
        ----------
        filename : str
            Output filename.
        **kwargs
            Additional arguments (not used for PyQtGraph).
        """
        if self._window is None:
            raise ValueError("Window must be created before saving")

        # PyQtGraph doesn't have a direct save method, so we use export
        exporter = pg.exporters.ImageExporter(self._plot)
        exporter.export(filename)

    @staticmethod
    def _color_to_hex(color: str) -> str:
        """
        Convert color name to hex format.

        Parameters
        ----------
        color : str
            Color name or hex string.

        Returns
        -------
        str
            Hex color string.
        """
        # If already hex, return as is
        if color.startswith('#'):
            return color

        # Common color mappings
        color_map = {
            'lightblue': '#ADD8E6',
            'lightgreen': '#90EE90',
            'salmon': '#FA8072',
            'gray': '#808080',
            'grey': '#808080',
            'red': '#FF0000',
            'blue': '#0000FF',
            'green': '#008000',
            'black': '#000000',
            'white': '#FFFFFF',
        }

        return color_map.get(color.lower(), '#000000')


# Register the backend
if HAS_PYQTGAPH:
    register_backend('pyqtgraph', PyQtGraphBackend)

