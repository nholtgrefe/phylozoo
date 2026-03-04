"""
Graph features module.

This module provides functions to extract and identify features of DirectedMultiGraph
instances (e.g., connectivity, self-loops, etc.).
"""

from typing import Any, Iterator, TypeVar

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


def bi_edge_connected_components(graph: 'DirectedMultiGraph') -> Iterator[set[T]]:
    """
    Get bi-edge connected components (2-edge-connected components) of the 
    the underlying undirected graph.
    
    A bi-edge connected component is a maximal subgraph that remains connected
    after removing any single edge. This is equivalent to finding connected components
    after removing all bridges (cut edges).
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each bi-edge connected component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> # Graph with cycle: 1-2-3-1, edge 3-4 is bridge
    >>> # Bi-edge connected components: {1, 2, 3}, {4}
    >>> G = DirectedMultiGraph()
    >>> _ = G.add_edge(1, 2)
    >>> _ = G.add_edge(2, 3)
    >>> _ = G.add_edge(3, 1)
    >>> _ = G.add_edge(3, 4)
    >>> comps = list(bi_edge_connected_components(G))
    >>> {1, 2, 3} in comps
    True
    >>> {4} in comps
    True
    >>> # Graph with two cycles connected by bridge: 1-2-3-1 and 4-5-6-4, connected via bridge 3-4
    >>> # Bi-edge connected components: {1, 2, 3}, {4, 5, 6}
    >>> G2 = DirectedMultiGraph()
    >>> _ = G2.add_edge(1, 2)
    >>> _ = G2.add_edge(2, 3)
    >>> _ = G2.add_edge(3, 1)
    >>> _ = G2.add_edge(3, 4)  # Bridge
    >>> _ = G2.add_edge(4, 5)
    >>> _ = G2.add_edge(5, 6)
    >>> _ = G2.add_edge(6, 4)
    >>> comps2 = list(bi_edge_connected_components(G2))
    >>> {1, 2, 3} in comps2
    True
    >>> {4, 5, 6} in comps2
    True
    """
    # Find all bridges (cut edges)
    bridges = set(nx.bridges(graph._combined))
    
    # Create a copy of the graph without bridges
    graph_without_bridges = graph._combined.copy()
    for u, v in bridges:
        # Remove all edges between u and v (handles parallel edges)
        # Note: bridges can't have parallel edges, but we remove all just to be safe
        while graph_without_bridges.has_edge(u, v):
            graph_without_bridges.remove_edge(u, v)
    
    # Return connected components of the graph without bridges
    return nx.connected_components(graph_without_bridges)


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


def cut_edges(
    graph: 'DirectedMultiGraph', 
    keys: bool = False, 
    data: bool | str = False
) -> set[tuple[T, T]] | set[tuple[T, T, int]] | set[tuple[T, T, Any]] | set[tuple[T, T, int, Any]] | list[tuple[T, T, dict[str, Any]]] | list[tuple[T, T, int, dict[str, Any]]]:
    """
    Find all cut-edges (bridges) in the graph.
    
    A cut-edge is an edge whose removal increases the number of
    weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    keys : bool, optional
        If True, return 3-tuples (u, v, key). If False, return 2-tuples (u, v).
        Default is False.
    data : bool | str, optional
        If False, return edges without data. If True, return edges with full data dict.
        If a string, return edges with the value of that attribute. Default is False.
    
    Returns
    -------
    set or list
        Cut-edges. Format depends on keys and data parameters:

        - keys=False, data=False: {(u, v), ...} (set)
        - keys=True, data=False: {(u, v, key), ...} (set)
        - keys=False, data=True: [(u, v, data_dict), ...] (list, since dicts are unhashable)
        - keys=True, data=True: [(u, v, key, data_dict), ...] (list, since dicts are unhashable)
        - keys=False, data='attr': {(u, v, attr_value), ...} (set)
        - keys=True, data='attr': {(u, v, key, attr_value), ...} (set)
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.features import cut_edges
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> edges = cut_edges(G)
    >>> (1, 2) in edges and (2, 3) in edges
    True
    >>> edges_with_keys = cut_edges(G, keys=True)
    >>> (1, 2, 0) in edges_with_keys and (2, 3, 0) in edges_with_keys
    True
    
    Notes
    -----
    This function uses Tarjan's algorithm for finding bridges, which runs in O(V + E) time.
    The algorithm works on the underlying undirected (weakly connected) representation.
    Parallel edges are never bridges, so this implementation optimizes by iterating through
    edges once and checking bridge membership.
    """
    # Get bridges from combined graph (O(V+E))
    bridges_set = set(nx.bridges(graph._combined))
    
    # Normalize bridges to (min, max) for efficient lookup
    bridges_normalized = {(min(u, v), max(u, v)) for u, v in bridges_set}
    
    # Use list for results with dicts (unhashable), set otherwise
    use_list = data is True
    result = [] if use_list else set()
    
    # Iterate through edges once and check if they're bridges
    # Bridges can't have parallel edges, so we only need to check each edge once
    for u, v, key, edge_data in graph._graph.edges(keys=True, data=True):
        # Check if this edge corresponds to a bridge
        edge_normalized = (min(u, v), max(u, v))
        if edge_normalized in bridges_normalized:
            # Format according to keys and data parameters
            if keys and data is True:
                result.append((u, v, key, edge_data))
            elif keys and isinstance(data, str):
                attr_val = edge_data.get(data)
                result.add((u, v, key, attr_val))
            elif keys:  # keys=True, data=False
                result.add((u, v, key))
            elif data is True:
                result.append((u, v, edge_data))
            elif isinstance(data, str):
                attr_val = edge_data.get(data)
                result.add((u, v, attr_val))
            else:  # keys=False, data=False
                result.add((u, v))
    
    return result


def cut_vertices(
    graph: 'DirectedMultiGraph',
    data: bool | str = False
) -> set[T] | set[tuple[T, Any]] | list[tuple[T, dict[str, Any]]]:
    """
    Find all cut-vertices (articulation points) in the graph.
    
    A cut-vertex is a vertex whose removal increases the number of
    weakly connected components.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to analyze.
    data : bool | str, optional
        If False, return vertices without data. If True, return vertices with full data dict.
        If a string, return vertices with the value of that attribute. Default is False.
    
    Returns
    -------
    set or list
        Cut-vertices. Format depends on data parameter:

        - data=False: {v, ...} (set)
        - data=True: [(v, data_dict), ...] (list, since dicts are unhashable)
        - data='attr': {(v, attr_value), ...} (set)
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.features import cut_vertices
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(2, 3)
    0
    >>> G.add_edge(2, 4)
    0
    >>> vertices = cut_vertices(G)
    >>> 2 in vertices
    True
    >>> 1 in vertices
    False
    
    Notes
    -----
    This function uses NetworkX's articulation_points algorithm, which runs in O(V + E) time.
    The algorithm works on the underlying undirected (weakly connected) representation.
    """
    # Get articulation points (O(V+E))
    art_points = set(nx.articulation_points(graph._combined))
    
    if data is False:
        return art_points
    
    # Use list for results with dicts (unhashable), set otherwise
    use_list = data is True
    result = [] if use_list else set()
    
    # Access node data directly from NetworkX graph (more efficient)
    nodes_data = graph._graph.nodes
    
    for v in art_points:
        if data is True:
            # Direct dict access is faster than creating a new dict
            node_data = dict(nodes_data[v]) if v in nodes_data else {}
            result.append((v, node_data))
        elif isinstance(data, str):
            # Direct attribute access
            node_data = nodes_data[v] if v in nodes_data else {}
            attr_val = node_data.get(data) if node_data else None
            result.add((v, attr_val))
    
    return result

