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
    """
    # Use NetworkX's bridge finding on the combined undirected graph
    # NetworkX's bridges() returns iterator of 2-tuples (u, v)
    bridges_2tuple = set(nx.bridges(graph._combined))
    
    # Use list for results with dicts (unhashable), set otherwise
    use_list = data is True
    result = [] if use_list else set()
    
    # For each bridge, find the corresponding directed edge(s)
    # Note: Parallel edges are never bridges, so we should only find single edges
    for u, v in bridges_2tuple:
        # Check both directions since _combined is undirected but _graph is directed
        edges_uv = []
        
        if u in graph._graph and v in graph._graph[u]:
            edge_dict = graph._graph[u][v]
            # Parallel edges are never bridges - skip if found (shouldn't happen)
            if len(edge_dict) > 1:
                continue  # This shouldn't happen; parallel edges can't be bridges
            for k, edge_data in edge_dict.items():
                edges_uv.append((u, v, k, edge_data))
        
        if v in graph._graph and u in graph._graph[v]:
            edge_dict = graph._graph[v][u]
            # Parallel edges are never bridges - skip if found (shouldn't happen)
            if len(edge_dict) > 1:
                continue  # This shouldn't happen; parallel edges can't be bridges
            for k, edge_data in edge_dict.items():
                edges_uv.append((v, u, k, edge_data))
        
        # Format according to keys and data parameters
        for u_dir, v_dir, k, edge_data in edges_uv:
            if keys and data is True:
                result.append((u_dir, v_dir, k, edge_data))
            elif keys and isinstance(data, str):
                attr_val = edge_data.get(data)
                result.add((u_dir, v_dir, k, attr_val))
            elif keys:  # keys=True, data=False
                result.add((u_dir, v_dir, k))
            elif data is True:
                result.append((u_dir, v_dir, edge_data))
            elif isinstance(data, str):
                attr_val = edge_data.get(data)
                result.add((u_dir, v_dir, attr_val))
            else:  # keys=False, data=False
                result.add((u_dir, v_dir))
    
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
    # Use NetworkX's articulation points on the combined undirected graph
    art_points = set(nx.articulation_points(graph._combined))
    
    if data is False:
        return art_points
    
    # Use list for results with dicts (unhashable), set otherwise
    use_list = data is True
    result = [] if use_list else set()
    
    for v in art_points:
        node_data = dict(graph._graph.nodes[v]) if v in graph._graph.nodes else {}
        if data is True:
            result.append((v, node_data))
        elif isinstance(data, str):
            attr_val = node_data.get(data) if node_data else None
            result.add((v, attr_val))
    
    return result

