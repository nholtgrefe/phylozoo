"""
Side module for directed level-k generators.

This module provides the Side class and its subclasses for representing
attachment points (sides) of generators.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

T = TypeVar('T')


@dataclass(frozen=True)
class Side:
    """
    Base class for sides (attachment points) of a generator.

    A side is a node or an edge where additional structure (leaves, trees) can be attached
    when building generators from lower-level generators.

    Subclasses: NodeSide (and HybridSide) for node attachment; EdgeSide and DirEdgeSide
    for edge attachment.
    """
    pass


# ---------------------------------------------------------------------------
# Node sides (attachment at a vertex)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NodeSide(Side):
    """
    Represents a node side of a generator (attachment at a vertex).

    Subclasses: IsolatedNodeSide for the single-vertex (level-0) generator;
    HybridSide for hybrid nodes (in-degree >= 2).

    Parameters
    ----------
    node : T
        The node identifier.

    Examples
    --------
    >>> node_side = NodeSide(node=0)
    >>> node_side.node
    0
    """
    node: T

    def __repr__(self) -> str:
        return f"NodeSide(node={self.node})"


@dataclass(frozen=True)
class IsolatedNodeSide(NodeSide):
    """
    Node side for the single vertex of a level-0 generator.

    A level-0 generator has one vertex and no edges. IsolatedNodeSide represents
    that vertex as the only attachment point (as opposed to HybridSide, which
    is for hybrid nodes). When attaching leaves for binary networks, exactly
    three leaves must be attached to this side.

    Parameters
    ----------
    node : T
        The node identifier (the unique vertex of the level-0 generator).

    Examples
    --------
    >>> side = IsolatedNodeSide(node=0)
    >>> side.node
    0
    >>> isinstance(side, NodeSide)
    True
    """

    def __repr__(self) -> str:
        return f"IsolatedNodeSide(node={self.node})"


@dataclass(frozen=True)
class HybridSide(NodeSide):
    """
    Represents a hybrid node side of a generator.

    A hybrid side is a node with in-degree >= 2 and out-degree 0 that serves as
    an attachment point. HybridSide is a NodeSide (attachment at a vertex).

    Parameters
    ----------
    node : T
        The node identifier for this hybrid side.

    Examples
    --------
    >>> hybrid_side = HybridSide(node=3)
    >>> hybrid_side.node
    3
    >>> isinstance(hybrid_side, NodeSide)
    True
    """

    def __repr__(self) -> str:
        return f"HybridSide(node={self.node})"


# ---------------------------------------------------------------------------
# Edge sides (attachment along an edge)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EdgeSide(Side):
    """
    Base class for edge sides (attachment along an edge).

    DirEdgeSide is the concrete directed-edge implementation.
    """
    pass


@dataclass(frozen=True)
class DirEdgeSide(EdgeSide):
    """
    Represents a directed edge side of a generator.

    A directed edge side is an edge (possibly parallel) that serves as
    an attachment point. The edge is identified by its endpoints and key.

    Parameters
    ----------
    u : T
        The source node of the edge.
    v : T
        The target node of the edge.
    key : int
        The edge key for parallel edges.

    Examples
    --------
    >>> edge_side = DirEdgeSide(u=8, v=4, key=0)
    >>> edge_side.u
    8
    >>> edge_side.v
    4
    >>> edge_side.key
    0
    """
    u: T
    v: T
    key: int

    def __repr__(self) -> str:
        return f"DirEdgeSide(u={self.u}, v={self.v}, key={self.key})"
