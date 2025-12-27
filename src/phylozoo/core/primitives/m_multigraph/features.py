"""
Graph features module.

This module provides functions to extract and identify features of MixedMultiGraph
instances (e.g., connectivity, self-loops, source components, etc.).
"""

from typing import Any, Iterator, TypeVar

import networkx as nx

from . import MixedMultiGraph

T = TypeVar('T')


def number_of_connected_components(graph: 'MixedMultiGraph') -> int:
    """
    Return the number of weakly connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    int
        Number of connected components.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_directed_edge(3, 4)
    0
    >>> number_of_connected_components(G)
    2
    """
    return nx.number_connected_components(graph._combined)


def is_connected(graph: 'MixedMultiGraph') -> bool:
    """
    Check if graph is weakly connected.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to check.
    
    Returns
    -------
    bool
        True if graph is connected, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_directed_edge(2, 3)
    0
    >>> is_connected(G)
    True
    """
    return nx.is_connected(graph._combined)


def has_parallel_edges(graph: 'MixedMultiGraph') -> bool:
    """
    Check if the graph has any parallel edges.
    
    Parallel edges are multiple edges between the same pair of nodes,
    either as multiple directed edges in the same direction or multiple undirected edges.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to check.
    
    Returns
    -------
    bool
        True if the graph has at least one pair of parallel edges, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> has_parallel_edges(G)
    False
    >>> G.add_directed_edge(1, 2)  # Add parallel directed edge
    1
    >>> has_parallel_edges(G)
    True
    """
    # Check directed edges
    for u, v in graph._directed.edges():
        if graph._directed.number_of_edges(u, v) > 1:
            return True
    
    # Check undirected edges
    for u, v in graph._undirected.edges():
        if graph._undirected.number_of_edges(u, v) > 1:
            return True
    
    return False


def connected_components(graph: 'MixedMultiGraph') -> Iterator[set[T]]:
    """
    Get weakly connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_directed_edge(3, 4)
    0
    >>> list(connected_components(G))
    [{1, 2}, {3, 4}]
    """
    return nx.connected_components(graph._combined)


def biconnected_components(graph: 'MixedMultiGraph') -> Iterator[set[T]]:
    """
    Get biconnected components of the underlying undirected graph.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each biconnected component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> # Single connected component with two biconnected components:
    >>> # Cycle 1: 1-2-3-1 and Cycle 2: 3-4-5-3 (articulation at 3)
    >>> _ = G.add_undirected_edge(1, 2)
    >>> _ = G.add_undirected_edge(2, 3)
    >>> _ = G.add_undirected_edge(3, 1)
    >>> _ = G.add_undirected_edge(3, 4)
    >>> _ = G.add_undirected_edge(4, 5)
    >>> _ = G.add_undirected_edge(5, 3)
    >>> comps = list(biconnected_components(G))
    >>> {1, 2, 3} in comps  # first cycle
    True
    >>> {3, 4, 5} in comps  # second cycle sharing articulation 3
    True
    """
    return nx.biconnected_components(graph._combined)


def bi_edge_connected_components(graph: 'MixedMultiGraph') -> Iterator[set[T]]:
    """
    Get bi-edge connected components (2-edge-connected components) of the graph.
    
    A bi-edge connected component is a maximal subgraph that remains connected
    after removing any single edge. This is equivalent to finding connected components
    after removing all bridges (cut edges).
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[set[T]]
        Iterator over sets of nodes in each bi-edge connected component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> # Graph with cycle: 1-2-3-1, edge 3-4 is bridge
    >>> # Bi-edge connected components: {1, 2, 3}, {4}
    >>> G = MixedMultiGraph()
    >>> _ = G.add_undirected_edge(1, 2)
    >>> _ = G.add_undirected_edge(2, 3)
    >>> _ = G.add_undirected_edge(3, 1)
    >>> _ = G.add_undirected_edge(3, 4)
    >>> comps = list(bi_edge_connected_components(G))
    >>> {1, 2, 3} in comps
    True
    >>> {4} in comps
    True
    >>> # Graph with two cycles connected by bridge: 1-2-3-1 and 4-5-6-4, connected via bridge 3-4
    >>> # Bi-edge connected components: {1, 2, 3}, {4, 5, 6}
    >>> G2 = MixedMultiGraph()
    >>> _ = G2.add_undirected_edge(1, 2)
    >>> _ = G2.add_undirected_edge(2, 3)
    >>> _ = G2.add_undirected_edge(3, 1)
    >>> _ = G2.add_undirected_edge(3, 4)  # Bridge
    >>> _ = G2.add_undirected_edge(4, 5)
    >>> _ = G2.add_undirected_edge(5, 6)
    >>> _ = G2.add_undirected_edge(6, 4)
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


def has_self_loops(graph: 'MixedMultiGraph') -> bool:
    """
    Check whether the mixed multigraph contains any self-loops (directed or undirected).
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to inspect.
    
    Returns
    -------
    bool
        True if at least one self-loop exists, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.features import has_self_loops
    >>> G = MixedMultiGraph()
    >>> has_self_loops(G)
    False
    >>> G.add_directed_edge(1, 1)
    0
    >>> has_self_loops(G)
    True
    >>> G.add_undirected_edge(2, 2)
    0
    >>> has_self_loops(G)
    True
    
    Notes
    -----
    MixedMultiGraph enforces mutual exclusivity between directed and undirected edges
    on the same node pair. Adding a directed self-loop to nodes that currently have an
    undirected self-loop will remove the undirected edge (and vice versa). This function
    checks both directed and undirected graphs, so it returns True if either side
    currently contains a self-loop.
    """
    return (
        nx.number_of_selfloops(graph._directed) > 0
        or nx.number_of_selfloops(graph._undirected) > 0
    )


def source_components(graph: 'MixedMultiGraph') -> list[tuple[list[T], list[tuple[T, T, int]], list[tuple[T, T, int]]]]:
    """
    Find all source components of a mixed multigraph.
    
    A source component is a connected component C of the undirected graph
    (i.e., the undirected part of the graph) with the property that there are no
    directed edges (u, v) in the multigraph with u not in C and v in C (i.e., no
    directed edges pointing into C).
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    list[tuple[list[T], list[tuple[T, T, int]], list[tuple[T, T, int]]]]
        For each source component, returns a tuple containing:
        - List of nodes in the component
        - List of undirected edges (u, v, key) within the component (includes all parallel edges)
        - List of directed edges (u, v, key) with u in the component and v not in the component
          (all outgoing edges of the component, includes all parallel edges)
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> G.add_directed_edge(3, 4)
    0
    >>> components = source_components(G)
    >>> len(components)
    1
    >>> nodes, undirected_edges, outgoing_edges = components[0]
    >>> sorted(nodes)
    [1, 2, 3]
    >>> sorted(undirected_edges)
    [(1, 2, 0), (2, 3, 0)]
    >>> outgoing_edges
    [(3, 4, 0)]
    """
    source_comps: list[tuple[list[T], list[tuple[T, T, int]], list[tuple[T, T, int]]]] = []
    
    # Get all connected components of the undirected graph
    undirected_components = list(nx.connected_components(graph._undirected))
    
    for component in undirected_components:
        component_set = set(component)
        
        # Check if there are any directed edges pointing into this component
        # (i.e., directed edges (u, v) where u is not in component and v is in component)
        has_incoming_edges = False
        for u, v, key in graph._directed.edges(keys=True):
            if u not in component_set and v in component_set:
                has_incoming_edges = True
                break
        
        # If no incoming edges, this is a source component
        if not has_incoming_edges:
            # Collect nodes in the component
            nodes = list(component)
            
            # Collect undirected edges within the component (including all parallel edges with keys)
            undirected_edges: list[tuple[T, T, int]] = []
            for u, v, key in graph._undirected.edges(keys=True):
                if u in component_set and v in component_set:
                    # Include all parallel edges (each with its key)
                    # Normalize to (min(u,v), max(u,v), key) to avoid duplicates from undirected representation
                    edge = (min(u, v), max(u, v), key)
                    undirected_edges.append(edge)
            
            # Collect directed edges (u, v, key) with u in component and v not in component
            # (includes all parallel edges with keys)
            outgoing_edges: list[tuple[T, T, int]] = []
            for u, v, key in graph._directed.edges(keys=True):
                if u in component_set and v not in component_set:
                    # Include all parallel edges (each with its key)
                    outgoing_edges.append((u, v, key))
            
            source_comps.append((nodes, undirected_edges, outgoing_edges))
    
    return source_comps


def cut_edges(
    graph: 'MixedMultiGraph', 
    keys: bool = False, 
    data: bool | str = False
) -> set[tuple[T, T]] | set[tuple[T, T, int]] | set[tuple[T, T, Any]] | set[tuple[T, T, int, Any]] | list[tuple[T, T, dict[str, Any]]] | list[tuple[T, T, int, dict[str, Any]]]:
    """
    Find all cut-edges (bridges) in the graph.
    
    A cut-edge is an edge whose removal increases the number of
    connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
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
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.features import cut_edges
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> edges = cut_edges(G)
    >>> (1, 2) in edges or (2, 1) in edges
    True
    >>> edges_with_keys = cut_edges(G, keys=True)
    >>> (1, 2, 0) in edges_with_keys or (2, 1, 0) in edges_with_keys
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
    
    # Track processed normalized edges to avoid duplicates
    processed_edges = set()
    
    # Check undirected edges
    for u, v, key, edge_data in graph._undirected.edges(keys=True, data=True):
        edge_normalized = (min(u, v), max(u, v))
        if edge_normalized in bridges_normalized and edge_normalized not in processed_edges:
            processed_edges.add(edge_normalized)
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
    
    # Check directed edges
    for u, v, key, edge_data in graph._directed.edges(keys=True, data=True):
        edge_normalized = (min(u, v), max(u, v))
        if edge_normalized in bridges_normalized and edge_normalized not in processed_edges:
            processed_edges.add(edge_normalized)
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
    graph: 'MixedMultiGraph',
    data: bool | str = False
) -> set[T] | set[tuple[T, Any]] | list[tuple[T, dict[str, Any]]]:
    """
    Find all cut-vertices (articulation points) in the graph.
    
    A cut-vertex is a vertex whose removal increases the number of
    connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
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
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.features import cut_vertices
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> G.add_undirected_edge(2, 4)
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
    
    # Access node data directly from NetworkX graphs (more efficient)
    directed_nodes = graph._directed.nodes
    undirected_nodes = graph._undirected.nodes
    
    for v in art_points:
        # For MixedMultiGraph, node data can be in either _directed or _undirected
        if data is True:
            # Direct dict access is faster than creating a new dict
            if v in directed_nodes:
                node_data = dict(directed_nodes[v])
            elif v in undirected_nodes:
                node_data = dict(undirected_nodes[v])
            else:
                node_data = {}
            result.append((v, node_data))
        elif isinstance(data, str):
            # Direct attribute access
            if v in directed_nodes:
                node_data = directed_nodes[v]
            elif v in undirected_nodes:
                node_data = undirected_nodes[v]
            else:
                node_data = {}
            attr_val = node_data.get(data) if node_data else None
            result.add((v, attr_val))
    
    return result


def _is_updown_path(graph: 'MixedMultiGraph', path: list[T], x: T, y: T) -> bool:
    """
    Check if a path is an up-down path from x to y.
    
    An up-down path is one where no two edges are oriented towards each other.
    Equivalently, it is a path where the first k edges can be oriented towards x
    and the remaining l-k edges can be oriented towards y.
    
    This implementation uses a two-pass algorithm:
    1. Find the turning point - the first vertex u_i where we encounter a directed
       edge (u_i, u_{i+1}) oriented from u_i towards y (i.e., u_i -> u_{i+1})
    2. After the turning point, no edges can be oriented towards x (i.e., no
       directed edge v->u where we traverse from u to v after the turning point)
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph.
    path : list[T]
        List of vertices forming a path from x to y. Assumed to be a valid path
        (i.e., edges exist between consecutive vertices).
    x : T
        Source vertex (should be path[0]).
    y : T
        Target vertex (should be path[-1]).
    
    Returns
    -------
    bool
        True if the path is an up-down path, False otherwise.
    
    Notes
    -----
    The input path is assumed to be valid (edges exist between consecutive vertices).
    """
    if len(path) < 2:
        return True
    
    turning_point_idx: int | None = None
    
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        if graph._directed.has_edge(u, v):
            turning_point_idx = i
            break
    
    if turning_point_idx is None:
        return True
    
    for i in range(turning_point_idx + 1, len(path) - 1):
        u, v = path[i], path[i + 1]
        if graph._directed.has_edge(v, u):
            return False
    
    return True


def updown_path_vertices(graph: 'MixedMultiGraph', x: T, y: T) -> set[T]:
    """
    Find all vertices on up-down paths between two vertices x and y.
    
    An up-down path from x to y is a path where no two edges are oriented towards
    each other. Equivalently, it is a path where the first k edges can be oriented
    towards x (where undirected edges can be oriented in either way) and the
    remaining l-k edges can be oriented towards y.
    
    This function finds all vertices v such that there exists an up-down path
    from x to y that passes through v.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph to analyze.
    x : T
        Source vertex.
    y : T
        Target vertex.
    
    Returns
    -------
    set[T]
        Set of all vertices on up-down paths from x to y, including x and y.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.features import updown_path_vertices
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> G.add_directed_edge(3, 4)
    0
    >>> G.add_undirected_edge(4, 5)
    0
    >>> vertices = updown_path_vertices(G, 1, 5)
    >>> vertices == {1, 2, 3, 4, 5}
    True
    
    Notes
    -----
    This implementation uses NetworkX's all_simple_paths to find all paths between
    x and y, then filters for up-down paths.
    """
    import networkx as nx
    
    if x not in graph.nodes() or y not in graph.nodes():
        return set()
    
    if x == y:
        return {x}
    
    # Use NetworkX to find all simple paths in the combined graph
    # The combined graph treats all edges as undirected for path finding
    # Use generator to avoid storing all paths in memory at once
    try:
        all_paths = nx.all_simple_paths(graph._combined, x, y)
    except nx.NetworkXNoPath:
        return set()
    
    # Collect all vertices on up-down paths
    result: set[T] = set()
    
    for path in all_paths:
        if _is_updown_path(graph, path, x, y):
            result.update(path)
            # Early termination: if we've found all vertices in the graph, we're done
            # (This is unlikely but can help in some cases)
            if len(result) == len(graph.nodes()):
                break
    
    return result

