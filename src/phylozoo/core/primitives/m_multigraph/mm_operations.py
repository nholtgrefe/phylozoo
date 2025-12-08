"""
Operations for MixedMultiGraph.

This module provides functions for working with MixedMultiGraph instances,
following NetworkX-style function-based API.
"""

from typing import TYPE_CHECKING, Iterator, List, Set, Tuple, TypeVar

import networkx as nx

if TYPE_CHECKING:
    from ..d_multigraph.dm_graph import DirectedMultiGraph
    from .mm_graph import MixedMultiGraph
else:
    from .mm_graph import MixedMultiGraph

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
    >>> from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph
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
    >>> from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_directed_edge(2, 3)
    0
    >>> is_connected(G)
    True
    """
    return nx.is_connected(graph._combined)


def connected_components(graph: 'MixedMultiGraph') -> Iterator[Set[T]]:
    """
    Get weakly connected components.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to analyze.
    
    Returns
    -------
    Iterator[Set[T]]
        Iterator over sets of nodes in each component.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_directed_edge(3, 4)
    0
    >>> list(connected_components(G))
    [{1, 2}, {3, 4}]
    """
    return nx.connected_components(graph._combined)


def identify_two_nodes(graph: 'MixedMultiGraph', u: T, v: T) -> None:
    """
    Identify two nodes u and v by keeping node u.
    
    Modifies the graph in place. All edges incident to v are moved to u,
    and v is removed.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to modify.
    u : T
        Node to keep.
    v : T
        Node to identify with u (will be removed).
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> identify_two_nodes(G, 1, 2)
    >>> list(G.nodes())
    [1, 3]
    """
    # Use NetworkX's contracted_nodes on each graph
    nx.contracted_nodes(graph._undirected, u, v, self_loops=False, copy=False)
    nx.contracted_nodes(graph._directed, u, v, self_loops=False, copy=False)
    nx.contracted_nodes(graph._combined, u, v, self_loops=False, copy=False)
    
    # Clean up any self-loops that might have been created
    # (contracted_nodes with self_loops=False should handle this, but be safe)
    if graph._undirected.has_edge(u, u):
        for k in list(graph._undirected[u][u].keys()):
            graph._undirected.remove_edge(u, u, key=k)
            graph._combined.remove_edge(u, u, key=k)
    if graph._directed.has_edge(u, u):
        for k in list(graph._directed[u][u].keys()):
            graph._directed.remove_edge(u, u, key=k)
            graph._combined.remove_edge(u, u, key=k)


def identify_node_set(graph: 'MixedMultiGraph', nodes: List[T] | Set[T]) -> None:
    """
    Identify all nodes in the set by keeping the first node.
    
    Modifies the graph in place.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to modify.
    nodes : List[T] | Set[T]
        Iterable of nodes to identify. The first node will be kept.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> G.add_undirected_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> identify_node_set(G, [1, 2, 3])
    >>> len(G.nodes()) <= 3
    True
    """
    nodes_list = list(nodes)
    if len(nodes_list) < 2:
        return
    
    for i in range(1, len(nodes_list)):
        identify_two_nodes(graph, nodes_list[0], nodes_list[i])


def source_components(graph: 'MixedMultiGraph') -> List[Tuple[List[T], List[Tuple[T, T, int]], List[Tuple[T, T, int]]]]:
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
    List[Tuple[List[T], List[Tuple[T, T, int]], List[Tuple[T, T, int]]]]
        For each source component, returns a tuple containing:
        - List of nodes in the component
        - List of undirected edges (u, v, key) within the component (includes all parallel edges)
        - List of directed edges (u, v, key) with u in the component and v not in the component
          (all outgoing edges of the component, includes all parallel edges)
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph
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
    source_comps: List[Tuple[List[T], List[Tuple[T, T, int]], List[Tuple[T, T, int]]]] = []
    
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
            undirected_edges: List[Tuple[T, T, int]] = []
            for u, v, key in graph._undirected.edges(keys=True):
                if u in component_set and v in component_set:
                    # Include all parallel edges (each with its key)
                    # Normalize to (min(u,v), max(u,v), key) to avoid duplicates from undirected representation
                    edge = (min(u, v), max(u, v), key)
                    undirected_edges.append(edge)
            
            # Collect directed edges (u, v, key) with u in component and v not in component
            # (includes all parallel edges with keys)
            outgoing_edges: List[Tuple[T, T, int]] = []
            for u, v, key in graph._directed.edges(keys=True):
                if u in component_set and v not in component_set:
                    # Include all parallel edges (each with its key)
                    outgoing_edges.append((u, v, key))
            
            source_comps.append((nodes, undirected_edges, outgoing_edges))
    
    return source_comps


