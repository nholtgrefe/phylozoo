"""
Base layout classes for SemiDirectedPhyNetwork.

This module provides layout data classes that store computed layout information
in a backend-agnostic format.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeVar

from ...utils.types import EdgeRoute, EdgeRoutes, Positions

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

T = TypeVar('T')


@dataclass(frozen=True)
class RadialLayout:
    """
    Radial layout for SemiDirectedPhyNetwork trees.

    This layout class stores the computed positions and routing information
    for a radial (circular) layout. The layout positions nodes in a circular
    arrangement with the root at the center and leaves on the outer circle.

    Attributes
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network this layout is for (read-only).
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
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork.layout import compute_radial_layout
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_radial_layout(net)
    >>> layout.get_position(1)
    (1.0, 0.0)
    >>> layout.algorithm
    'radial'
    """

    network: 'SemiDirectedPhyNetwork'
    positions: Positions[T]
    edge_routes: EdgeRoutes[T]
    algorithm: str = 'radial'
    parameters: dict[str, Any] = field(default_factory=dict)

    def get_position(self, node: T) -> tuple[float, float]:
        """
        Get position of a node.

        Parameters
        ----------
        node : T
            Node ID.

        Returns
        -------
        tuple[float, float]
            (x, y) position.

        Raises
        ------
        KeyError
            If node is not in the layout.

        Examples
        --------
        >>> layout = compute_radial_layout(network)
        >>> x, y = layout.get_position(1)
        >>> isinstance(x, float)
        True
        """
        return self.positions[node]

    def get_edge_route(
        self, u: T, v: T, key: int = 0
    ) -> 'EdgeRoute | None':
        """
        Get routing information for an edge.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int, optional
            Edge key (for parallel edges). By default 0.

        Returns
        -------
        EdgeRoute | None
            Edge routing information, or None if edge not found.

        Examples
        --------
        >>> layout = compute_radial_layout(network)
        >>> route = layout.get_edge_route(1, 2, 0)
        >>> route is not None
        True
        """
        return self.edge_routes.get((u, v, key))

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"RadialLayout(network={self.network}, "
            f"nodes={len(self.positions)}, edges={len(self.edge_routes)}, "
            f"algorithm={self.algorithm})"
        )
