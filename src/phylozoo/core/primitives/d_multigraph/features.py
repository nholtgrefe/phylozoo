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


def has_parallel_edges(graph: 'DirectedMultiGraph') -> bool:
    """
    Check if the graph has any parallel edges.
    
    Parallel edges are multiple edges between the same pair of nodes in the same direction.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to check.
    
    Returns
    -------
    bool
        True if the graph has at least one pair of parallel edges, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> has_parallel_edges(G)
    False
    >>> G.add_edge(1, 2)  # Add parallel edge
    1
    >>> has_parallel_edges(G)
    True
    """
    for u, v in graph._graph.edges():
        if graph._graph.number_of_edges(u, v) > 1:
            return True
    return False


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


def biconnected_components(graph: 'DirectedMultiGraph') -> Iterator[set[T]]:
    """
    Get biconnected components of the underlying undirected graph.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each biconnected component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> # Single connected component with two biconnected components:
    >>> # Cycle 1: 1-2-3-1 and Cycle 2: 3-4-5-3 (articulation at 3)
    >>> _ = G.add_edge(1, 2)
    >>> _ = G.add_edge(2, 3)
    >>> _ = G.add_edge(3, 1)
    >>> _ = G.add_edge(3, 4)
    >>> _ = G.add_edge(4, 5)
    >>> _ = G.add_edge(5, 3)
    >>> comps = list(biconnected_components(G))
    >>> {1, 2, 3} in comps  # first cycle
    True
    >>> {3, 4, 5} in comps  # second cycle sharing articulation 3
    True
    """
    return nx.biconnected_components(graph._combined)


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


def is_cutedge(graph: 'DirectedMultiGraph', u: T, v: T, key: int | None = None) -> bool:
    """
    Check if edge (u, v) is a cut-edge.
    
    A cut-edge is an edge whose removal increases the number of
    weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    u : T
        Source node.
    v : T
        Target node.
    key : int | None, optional
        Edge key. If None, checks the first edge found. By default None.
    
    Returns
    -------
    bool
        True if (u, v) is a cut-edge, False otherwise.
    
    Raises
    ------
    ValueError
        If (u, v) is not an edge in the graph.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.features import is_cutedge
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> is_cutedge(G, 2, 3)
    True
    """
    if not graph.has_edge(u, v, key):
        raise ValueError(f"Edge ({u}, {v}, {key}) is not in the graph.")
    
    # Find edge key if not provided
    if key is None:
        if u in graph._graph and v in graph._graph[u]:
            key = next(iter(graph._graph[u][v].keys()))
        else:
            raise ValueError(f"Edge ({u}, {v}) is not in the graph.")
    
    num_before = nx.number_connected_components(graph._combined)
    
    # Temporarily remove and check
    graph._combined.remove_edge(u, v, key)
    num_after = nx.number_connected_components(graph._combined)
    
    # Restore edge
    edge_data = dict(graph._graph[u][v][key])
    graph._combined.add_edge(u, v, key=key, **edge_data)
    
    return num_after != num_before


def is_cutvertex(graph: 'DirectedMultiGraph', v: T) -> bool:
    """
    Check if v is a cut-vertex.
    
    A cut-vertex is a vertex whose removal increases the number of
    weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    v : T
        Vertex.
    
    Returns
    -------
    bool
        True if v is a cut-vertex, False otherwise.
    
    Raises
    ------
    ValueError
        If v is not a vertex in the graph.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.features import is_cutvertex
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> is_cutvertex(G, 2)
    True
    """
    if v not in graph:
        raise ValueError(f"Vertex {v} is not in the graph.")
    
    num_before = nx.number_connected_components(graph._combined)
    
    # Temporarily remove and check
    graph._combined.remove_node(v)
    num_after = nx.number_connected_components(graph._combined)
    
    # Restore node and its edges
    graph._combined.add_node(v)
    # Restore edges incident to v
    for predecessor in graph._graph.predecessors(v):
        for k, data in graph._graph[predecessor][v].items():
            graph._combined.add_edge(predecessor, v, key=k, **data)
    for successor in graph._graph.successors(v):
        for k, data in graph._graph[v][successor].items():
            graph._combined.add_edge(v, successor, key=k, **data)
    
    return num_after != num_before

