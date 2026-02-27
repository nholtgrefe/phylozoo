"""
Side module for semi-directed level-k generators.

This module provides UndirEdgeSide (subclass of EdgeSide) for undirected edge sides,
BidirectedEdgeSide (subclass of EdgeSide) for the level-1 bidirected self-loop,
and re-exports Side, NodeSide, IsolatedNodeSide, HybridSide, EdgeSide, and DirEdgeSide
from the directed generator module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

from ...dnetwork.generator.side import (
    DirEdgeSide,
    EdgeSide,
    HybridSide,
    IsolatedNodeSide,
    NodeSide,
    Side,
)

T = TypeVar("T")


@dataclass(frozen=True)
class UndirEdgeSide(EdgeSide):
    """
    Represents an undirected edge side of a generator.

    An undirected edge side represents an undirected edge (possibly parallel) that
    serves as an attachment point. The edge is identified by its endpoints and key.

    Parameters
    ----------
    u : T
        The first node of the edge.
    v : T
        The second node of the edge.
    key : int
        The edge key for parallel edges.

    Examples
    --------
    >>> edge_side = UndirEdgeSide(u=1, v=2, key=0)
    >>> edge_side.u
    1
    >>> edge_side.v
    2
    >>> edge_side.key
    0
    """
    u: T
    v: T
    key: int

    def __repr__(self) -> str:
        return f"UndirEdgeSide(u={self.u}, v={self.v}, key={self.key})"


@dataclass(frozen=True)
class BidirectedEdgeSide(EdgeSide):
    """
    Represents the bidirected edge of a level-1 semi-directed generator.

    A level-1 semi-directed generator is a single node with an undirected self-loop
    (bidirected edge). This cannot be represented as a directed self-loop in the
    underlying graph, so it is stored as an undirected self-loop. BidirectedEdgeSide
    identifies that edge by the node and the edge key.

    Parameters
    ----------
    node : T
        The single node incident to the undirected self-loop.
    key : int
        The edge key for the undirected self-loop.

    Examples
    --------
    >>> side = BidirectedEdgeSide(node=0, key=0)
    >>> side.node
    0
    >>> side.key
    0
    """
    node: T
    key: int

    def __repr__(self) -> str:
        return f"BidirectedEdgeSide(node={self.node}, key={self.key})"
