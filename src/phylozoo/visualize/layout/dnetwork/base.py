"""
Base layout classes for DirectedPhyNetwork.

This module provides layout classes specifically for DirectedPhyNetwork.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeVar

from ..base import Layout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


@dataclass(frozen=True)
class DNetLayout(Layout):
    """
    Layout class for DirectedPhyNetwork.

    Extends the base Layout class to store a DirectedPhyNetwork instead of
    a generic graph object. This provides network-specific functionality
    while maintaining the general layout interface.

    Attributes
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network this layout is for (read-only).
    positions : dict[T, tuple[float, float]]
        Node positions mapping node ID to (x, y) coordinates (read-only).
    edge_routes : dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute (read-only).
    algorithm : str
        Name of the layout algorithm used (read-only).
    parameters : dict[str, Any]
        Parameters used to generate this layout (read-only).

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.layout.dnetwork import DNetLayout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> # Layouts are typically created via compute functions
    >>> # layout = compute_rectangular_dnet_layout(net)
    """

    # network is stored as graph in parent, expose it as network property
    @property
    def network(self) -> 'DirectedPhyNetwork':
        """Get the network (stored as graph in parent)."""
        return self.graph

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"DNetLayout(network={self.network}, "
            f"nodes={len(self.positions)}, edges={len(self.edge_routes)}, "
            f"algorithm={self.algorithm})"
        )


@dataclass(frozen=True)
class RectangularDNetLayout(DNetLayout):
    """
    Rectangular layout for DirectedPhyNetwork.

    Specialized layout class for rectangular (layered) layouts of
    DirectedPhyNetwork. Stores additional information specific to
    rectangular layouts such as layer assignments and orientation.

    Attributes
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network this layout is for (read-only).
    positions : dict[T, tuple[float, float]]
        Node positions mapping node ID to (x, y) coordinates (read-only).
    edge_routes : dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute (read-only).
    layers : dict[T, int]
        Layer assignment for each node (0 = root layer) (read-only).
    orientation : str
        Layout orientation: 'top-bottom' or 'left-right' (read-only).
    algorithm : str
        Name of the layout algorithm used (read-only).
    parameters : dict[str, Any]
        Parameters used to generate this layout (read-only).

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.layout.dnetwork import compute_rectangular_dnet_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_rectangular_dnet_layout(net)
    >>> layout.get_position(1)
    (0.0, 1.5)
    >>> layout.get_layer(1)
    1
    """

    layers: dict[T, int] = field(default_factory=dict)
    orientation: str = 'top-bottom'

    def get_layer(self, node: T) -> int:
        """
        Get layer assignment for a node.

        Parameters
        ----------
        node : T
            Node ID.

        Returns
        -------
        int
            Layer number (0 = root layer).

        Raises
        ------
        KeyError
            If node is not in the layout.

        Examples
        --------
        >>> layout = compute_rectangular_dnet_layout(net)
        >>> layer = layout.get_layer(1)
        >>> layer >= 0
        True
        """
        return self.layers[node]

    def to_dict(self) -> dict[str, Any]:
        """
        Convert layout to dictionary (for serialization).

        Returns
        -------
        dict[str, Any]
            Dictionary representation of the layout.
        """
        from typing import Any
        data = super().to_dict()
        data['layers'] = self.layers
        data['orientation'] = self.orientation
        return data

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"RectangularDNetLayout(network={self.network}, "
            f"nodes={len(self.positions)}, edges={len(self.edge_routes)}, "
            f"algorithm={self.algorithm}, orientation={self.orientation})"
        )

