"""
Side module for semi-directed level-k generators.

This module provides the UndirEdgeSide class for representing undirected edge sides,
and re-exports the base Side classes from the directed generator module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

from ...dnetwork.generator.side import DirEdgeSide, HybridSide, Side

T = TypeVar('T')


@dataclass(frozen=True)
class UndirEdgeSide(Side):
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

