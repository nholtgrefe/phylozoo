"""
Base rendering classes for DirectedPhyNetwork.

This module provides abstract base classes for rendering network layouts.
Actual rendering is performed by backend implementations.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..layout.base import DAGLayout
    from ..styling.style import NetworkStyle


class Renderer(ABC):
    """
    Abstract base class for rendering network layouts.

    This class defines the interface for rendering network layouts.
    Backend-specific implementations should inherit from this class.

    Notes
    -----
    In practice, rendering is handled by Backend implementations.
    This class is provided for potential future use or alternative
    rendering architectures.
    """

    @abstractmethod
    def render_edges(
        self, layout: 'DAGLayout', style: 'NetworkStyle'
    ) -> dict[str, Any]:
        """
        Render all edges in the layout.

        Parameters
        ----------
        layout : DAGLayout
            The layout to render.
        style : NetworkStyle
            Styling configuration.

        Returns
        -------
        dict[str, Any]
            Dictionary mapping edge keys to rendered elements.
        """
        pass

    @abstractmethod
    def render_nodes(
        self, layout: 'DAGLayout', style: 'NetworkStyle'
    ) -> dict[str, Any]:
        """
        Render all nodes in the layout.

        Parameters
        ----------
        layout : DAGLayout
            The layout to render.
        style : NetworkStyle
            Styling configuration.

        Returns
        -------
        dict[str, Any]
            Dictionary mapping node IDs to rendered elements.
        """
        pass

