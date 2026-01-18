"""
Base layout classes for DirectedMultiGraph.

This module provides layout data classes that store computed layout information.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, TypeVar

from phylozoo.viz._types import EdgeRoutes, Layout, Positions

if TYPE_CHECKING:
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph

T = TypeVar('T')


@dataclass(frozen=True)
class DMGraphLayout(Layout[T]):
    """
    Layout for DirectedMultiGraph.

    This layout class stores the computed positions and routing information
    for a DirectedMultiGraph layout.

    Attributes
    ----------
    network : DirectedMultiGraph
        The directed multigraph this layout is for (read-only).
    positions : dict[T, tuple[float, float]]
        Node positions mapping node ID to (x, y) coordinates (read-only).
    edge_routes : dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute (read-only).
    algorithm : str
        Name of the layout algorithm used (read-only).
    parameters : dict[str, Any]
        Parameters used to generate this layout (read-only).
    """

    network: 'DirectedMultiGraph[T]'
    positions: Positions[T]
    edge_routes: EdgeRoutes[T]
    algorithm: str = 'unknown'
    parameters: dict[str, Any] = field(default_factory=dict)
