"""
Graph features module.

This module provides functions to extract and identify features of DirectedMultiGraph
instances (e.g., connectivity, self-loops, etc.).
"""

from typing import Iterator, TypeVar

import networkx as nx

from . import DirectedMultiGraph

T = TypeVar('T')


def number_of_connected_components(graph: 'DirectedMultiGraph') -> int:
    """
    Return the number of weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    int
        Number of connected components.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(3, 4)
    0
    >>> number_of_connected_components(G)
    2
    """
    return nx.number_connected_components(graph._combined)


def is_connected(graph: 'DirectedMultiGraph') -> bool:
    """
    Check if graph is weakly connected.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to check.
    
    Returns
    -------
    bool
        True if graph is connected, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> is_connected(G)
    True
    """
    return nx.is_connected(graph._combined)


def connected_components(graph: 'DirectedMultiGraph') -> Iterator[set[T]]:
    """
    Get weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(3, 4)
    0
    >>> list(connected_components(G))
    [{1, 2}, {3, 4}]
    """
    return nx.connected_components(graph._combined)


def has_self_loops(graph: 'DirectedMultiGraph') -> bool:
    """
    Check whether the directed multigraph contains any self-loops.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to inspect.
    
    Returns
    -------
    bool
        True if at least one self-loop exists, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.features import has_self_loops
    >>> G = DirectedMultiGraph()
    >>> has_self_loops(G)
    False
    >>> G.add_edge(1, 1)
    0
    >>> has_self_loops(G)
    True
    """
    return nx.number_of_selfloops(graph._graph) > 0

