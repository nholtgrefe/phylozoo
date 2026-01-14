"""
Backend abstraction for DirectedPhyNetwork plotting.

This module provides the abstract base class for plotting backends and
a registry system for managing different backend implementations.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ...utils.types import EdgeRoute
    from ..layout.base import DAGLayout
    from ..styling.style import NetworkStyle

# Backend registry
_backends: dict[str, type['Backend']] = {}


def register_backend(name: str, backend_class: type['Backend']) -> None:
    """
    Register a backend implementation.

    Parameters
    ----------
    name : str
        Name of the backend (e.g., 'matplotlib', 'pyqtgraph').
    backend_class : type[Backend]
        Backend class that implements the Backend interface.

    Examples
    --------
    >>> from phylozoo.viz.dnetwork.backends import register_backend
    >>> register_backend('my_backend', MyBackend)
    """
    _backends[name] = backend_class


def get_backend(name: str) -> type['Backend']:
    """
    Get a backend class by name.

    Parameters
    ----------
    name : str
        Name of the backend.

    Returns
    -------
    type[Backend]
        Backend class.

    Raises
    ------
    ValueError
        If backend is not registered.

    Examples
    --------
    >>> from phylozoo.viz.dnetwork.backends import get_backend
    >>> BackendClass = get_backend('matplotlib')
    """
    if name not in _backends:
        raise ValueError(
            f"Backend '{name}' is not registered. "
            f"Available backends: {list(_backends.keys())}"
        )
    return _backends[name]


class Backend(ABC):
    """
    Abstract base class for plotting backends.

    This class defines the interface that all plotting backends must implement.
    Backends are responsible for rendering network layouts to specific output
    formats (matplotlib figures, pyqtgraph windows, etc.).

    Attributes
    ----------
    figure : Any
        Backend-specific figure object (read-only).
    axes : Any
        Backend-specific axes/plot area object (read-only).
    """

    def __init__(self) -> None:
        """Initialize the backend."""
        self._figure: Any = None
        self._axes: Any = None

    @property
    def figure(self) -> Any:
        """
        Get the figure object.

        Returns
        -------
        Any
            Backend-specific figure object.
        """
        return self._figure

    @property
    def axes(self) -> Any:
        """
        Get the axes object.

        Returns
        -------
        Any
            Backend-specific axes object.
        """
        return self._axes

    @abstractmethod
    def create_figure(self, **kwargs: Any) -> Any:
        """
        Create a new figure.

        Parameters
        ----------
        **kwargs
            Backend-specific figure creation parameters.

        Returns
        -------
        Any
            Backend-specific figure object.
        """
        pass

    @abstractmethod
    def create_axes(self, figure: Any | None = None, **kwargs: Any) -> Any:
        """
        Create axes on a figure.

        Parameters
        ----------
        figure : Any, optional
            Existing figure to create axes on. If None, creates a new figure.
        **kwargs
            Backend-specific axes creation parameters.

        Returns
        -------
        Any
            Backend-specific axes object.
        """
        pass

    @abstractmethod
    def render_edge(
        self,
        route: 'EdgeRoute',
        style: 'NetworkStyle',
        layout: 'DAGLayout',
        parallel_offset: float = 0.0,
    ) -> Any:
        """
        Render a single edge.

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
        Any
            Backend-specific rendered element.
        """
        pass

    @abstractmethod
    def render_node(
        self,
        node: Any,
        position: tuple[float, float],
        node_type: str,
        style: 'NetworkStyle',
    ) -> Any:
        """
        Render a single node.

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
        Any
            Backend-specific rendered element.
        """
        pass

    @abstractmethod
    def add_label(
        self,
        position: tuple[float, float],
        text: str,
        style: 'NetworkStyle',
        node_type: str = 'leaf',
    ) -> Any:
        """
        Add a text label.

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
        Any
            Backend-specific label element.
        """
        pass

    @abstractmethod
    def show(self) -> None:
        """
        Display the plot.

        This method should display the plot in an appropriate way for
        the backend (e.g., show matplotlib figure, display pyqtgraph window).
        """
        pass

    @abstractmethod
    def save(self, path: str, **kwargs: Any) -> None:
        """
        Save the plot to a file.

        Parameters
        ----------
        path : str
            File path to save to.
        **kwargs
            Backend-specific save parameters.
        """
        pass

