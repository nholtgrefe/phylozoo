"""
Base layout classes for SemiDirectedPhyNetwork.

This module provides layout data classes that store computed layout information.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeVar

from phylozoo.viz._types import Positions
from ...m_multigraph.layout.base import MGraphLayout

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

T = TypeVar('T')


@dataclass(frozen=True)
class SDNetLayout(MGraphLayout[T]):
    """
    Layout for SemiDirectedPhyNetwork.

    This layout class stores the computed positions and routing information
    for a SemiDirectedPhyNetwork layout. It is the base layout class for all
    sdnetwork layout algorithms.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork.layout import compute_pz_radial_layout
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_pz_radial_layout(net)
    >>> layout.get_position(1)
    (1.0, 0.0)
    >>> layout.algorithm
    'pz-radial'
    
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
    """

    network: 'SemiDirectedPhyNetwork'
    algorithm: str = 'pz-radial'  # Default, but typically overridden by specific layouts
