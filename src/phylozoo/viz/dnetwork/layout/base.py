"""
Base layout classes for DirectedPhyNetwork.

This module provides layout data classes that store computed layout information.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeVar

# Positions inherited from DMGraphLayout
from ...graphs.dmgraph.layout.base import DMGraphLayout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


@dataclass(frozen=True)
class DNetLayout(DMGraphLayout[T]):
    """
    Layout for DirectedPhyNetwork.

    This layout class stores the computed positions and routing information
    for a DirectedPhyNetwork layout. It is the base layout class for all
    dnetwork layout algorithms.

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
    >>> from phylozoo.viz.dnetwork.layout import compute_pz_dag_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_pz_dag_layout(net)
    >>> layout.get_position(1)
    (0.0, 1.5)
    >>> layout.algorithm
    'pz-dag'
    """

    network: 'DirectedPhyNetwork'
    backbone_edges: set[tuple[T, T, int]] = field(default_factory=set)
    reticulate_edges: set[tuple[T, T, int]] = field(default_factory=set)
    algorithm: str = 'pz-dag'  # Default, but typically overridden by specific layouts
