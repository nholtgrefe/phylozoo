"""
Graph features module.

This module provides functions to extract and identify features of MixedMultiGraph
instances (e.g., connectivity, self-loops, source components, etc.).
"""

from typing import Iterator, TypeVar

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


def is_cutedge(graph: 'MixedMultiGraph', u: T, v: T, key: int | None = None) -> bool:
    """
    Check if edge (u, v) is a cut-edge.
    
    A cut-edge is an edge whose removal increases the number of
    connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
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
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.features import is_cutedge
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> is_cutedge(G, 2, 3)
    True
    """
    if not graph.has_edge(u, v, key):
        raise ValueError(f"Edge ({u}, {v}, {key}) is not in the graph.")
    
    # Find edge key if not provided
    if key is None:
        if graph._undirected.has_edge(u, v):
            key = next(iter(graph._undirected[u][v].keys()))
        elif graph._directed.has_edge(u, v):
            key = next(iter(graph._directed[u][v].keys()))
        else:
            raise ValueError(f"Edge ({u}, {v}) is not in the graph.")
    
    num_before = nx.number_connected_components(graph._combined)
    
    # Temporarily remove and check
    graph._combined.remove_edge(u, v, key)
    num_after = nx.number_connected_components(graph._combined)
    
    # Restore edge
    if graph._undirected.has_edge(u, v, key):
        edge_data = dict(graph._undirected[u][v][key])
        graph._combined.add_edge(u, v, key=key, **edge_data)
    elif graph._directed.has_edge(u, v, key):
        edge_data = dict(graph._directed[u][v][key])
        graph._combined.add_edge(u, v, key=key, **edge_data)
    
    return num_after != num_before


def is_cutvertex(graph: 'MixedMultiGraph', v: T) -> bool:
    """
    Check if v is a cut-vertex.
    
    A cut-vertex is a vertex whose removal increases the number of
    connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
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
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.features import is_cutvertex
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
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
    if v in graph._undirected:
        for neighbor in graph._undirected.neighbors(v):
            for k, data in graph._undirected[v][neighbor].items():
                graph._combined.add_edge(v, neighbor, key=k, **data)
    if v in graph._directed:
        for predecessor in graph._directed.predecessors(v):
            for k, data in graph._directed[predecessor][v].items():
                graph._combined.add_edge(predecessor, v, key=k, **data)
        for successor in graph._directed.successors(v):
            for k, data in graph._directed[v][successor].items():
                graph._combined.add_edge(v, successor, key=k, **data)
    
    return num_after != num_before

