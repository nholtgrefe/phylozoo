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
    
    This is a general wrapper class. Use HybridSide for node sides and DirEdgeSide
    for edge sides.
    """
    pass


@dataclass(frozen=True)
class HybridSide(Side):
    """
    Represents a hybrid side of a generator.
    
    A hybrid side is a node with in-degree 2 and out-degree 0.
    This represents a hybrid node that serves as an attachment point.
    
    Parameters
    ----------
    node : T
        The node identifier for this hybrid side.
    
    Examples
    --------
    >>> hybrid_side = HybridSide(node=3)
    >>> hybrid_side.node
    3
    """
    node: T
    
    def __repr__(self) -> str:
        return f"HybridSide(node={self.node})"


@dataclass(frozen=True)
class DirEdgeSide(Side):
    """
    Represents a directed edge side of a generator.
    
    A directed edge side represents an edge (possibly parallel) that serves as
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

