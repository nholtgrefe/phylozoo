"""
Base layout classes for DirectedPhyNetwork.

This module provides layout data classes that store computed layout information
in a backend-agnostic format.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeVar

from ...utils.types import EdgeRoute, EdgeRoutes, Positions

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


@dataclass(frozen=True)
class DAGLayout:
    """
    DAG layout for DirectedPhyNetwork.

    This layout class stores the computed positions and routing information
    for a DAG layout. The layout uses a tree-backbone heuristic with crossing
    minimization.

    Attributes
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network this layout is for (read-only).
    positions : dict[T, tuple[float, float]]
        Node positions mapping node ID to (x, y) coordinates (read-only).
    edge_routes : dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute (read-only).
    backbone_edges : set[tuple[T, T, int]]
        Edges in the backbone tree (read-only).
    reticulate_edges : set[tuple[T, T, int]]
        Hybrid edges not in the backbone tree (read-only).
    algorithm : str
        Name of the layout algorithm used (read-only).
    parameters : dict[str, Any]
        Parameters used to generate this layout (read-only).

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz.dnetwork.layout import compute_dag_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_dag_layout(net)
    >>> layout.get_position(1)
    (0.0, 1.5)
    >>> layout.algorithm
    'dag'
    """

    network: 'DirectedPhyNetwork'
    positions: Positions[T]
    edge_routes: EdgeRoutes[T]
    backbone_edges: set[tuple[T, T, int]] = field(default_factory=set)
    reticulate_edges: set[tuple[T, T, int]] = field(default_factory=set)
    algorithm: str = 'dag'
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
        >>> layout = compute_dag_layout(network)
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
        >>> layout = compute_dag_layout(network)
        >>> route = layout.get_edge_route(1, 2, 0)
        >>> route is not None
        True
        """
        return self.edge_routes.get((u, v, key))

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"DAGLayout(network={self.network}, "
            f"nodes={len(self.positions)}, edges={len(self.edge_routes)}, "
            f"algorithm={self.algorithm})"
        )

