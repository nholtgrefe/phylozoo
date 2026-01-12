"""
Isomorphism module for mixed multi-graphs.

This module provides functions for checking graph isomorphism between
MixedMultiGraph instances.
"""

from __future__ import annotations

from typing import Any, TypeVar

import networkx as nx

from .base import MixedMultiGraph

T = TypeVar('T')


def _to_digraph_for_isomorphism(
    G: MixedMultiGraph[T],
) -> nx.MultiDiGraph:
    """
    Convert MixedMultiGraph to MultiDiGraph for isomorphism checking.
    
    Undirected edges are converted to bidirectional directed edges (u->v and v->u
    with the same key and attributes). This works because MixedMultiGraph enforces
    mutual exclusivity: edges between the same nodes are either all directed or
    all undirected, never mixed.
    
    Parameters
    ----------
    G : MixedMultiGraph
        Mixed multi-graph to convert.
    
    Returns
    -------
    nx.MultiDiGraph
        NetworkX MultiDiGraph with all edges as directed. Undirected edges appear
        as bidirectional pairs.
    """
    result = nx.MultiDiGraph()
    
    # Copy graph attributes
    if G._directed.graph:
        result.graph.update(G._directed.graph)
    
    # Copy nodes with attributes (from both subgraphs, merge if needed)
    for node, data in G._directed.nodes(data=True):
        result.add_node(node, **data)
    for node, data in G._undirected.nodes(data=True):
        if node not in result:
            result.add_node(node, **data)
        else:
            # Merge attributes if node exists in both
            result.nodes[node].update(data)
    
    # Add directed edges as-is
    for u, v, key, data in G._directed.edges(keys=True, data=True):
        result.add_edge(u, v, key=key, **data)
    
    # Add undirected edges as bidirectional directed edges
    # Use the same key for both directions to preserve parallel edge structure
    for u, v, key, data in G._undirected.edges(keys=True, data=True):
        result.add_edge(u, v, key=key, **data)
        result.add_edge(v, u, key=key, **data)
    
    return result


def is_isomorphic(
    G1: MixedMultiGraph[T],
    G2: MixedMultiGraph[T],
    node_attrs: list[str] | None = None,
    edge_attrs: list[str] | None = None,
    graph_attrs: list[str] | None = None,
) -> bool:
    """
    Check if two mixed multi-graphs are isomorphic.
    
    Two graphs are isomorphic if there exists a bijection between their node sets
    that preserves adjacency, edge direction, and parallel edges. Undirected edges
    are treated as bidirectional, so they must match with other undirected edges
    (which appear as bidirectional pairs in the converted graph).
    
    Parameters
    ----------
    G1 : MixedMultiGraph
        First mixed multi-graph.
    G2 : MixedMultiGraph
        Second mixed multi-graph.
    node_attrs : list[str] | None, optional
        List of node attribute names to match. If None, node attributes are ignored.
        Nodes must have matching values for all specified attributes to be considered
        isomorphic. If a node doesn't have an attribute, it matches only with nodes
        that also don't have that attribute. By default None.
    edge_attrs : list[str] | None, optional
        List of edge attribute names to match. If None, edge attributes are ignored.
        Edges must have matching values for all specified attributes. If a node 
        doesn't have an attribute, it matches only with nodes that also don't have 
        that attribute.By default None.
    graph_attrs : list[str] | None, optional
        List of graph-level attribute names to match. If None, graph attributes are
        ignored. Graph attributes are checked before isomorphism checking for efficiency.
        By default None.
    
    Returns
    -------
    bool
        True if the graphs are isomorphic, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> G1 = MixedMultiGraph(undirected_edges=[(1, 2)], directed_edges=[(2, 3)])
    >>> G2 = MixedMultiGraph(undirected_edges=[(4, 5)], directed_edges=[(5, 6)])
    >>> is_isomorphic(G1, G2)
    True
    >>> # With node labels
    >>> G1.add_node(1, label='A')
    >>> G2.add_node(4, label='A')
    >>> is_isomorphic(G1, G2, node_attrs=['label'])
    True
    >>> # Edge attributes
    >>> G3 = MixedMultiGraph(undirected_edges=[{'u': 1, 'v': 2, 'weight': 1.0}])
    >>> G4 = MixedMultiGraph(undirected_edges=[{'u': 3, 'v': 4, 'weight': 1.0}])
    >>> is_isomorphic(G3, G4, edge_attrs=['weight'])
    True
    >>> # Mixed directed and undirected edges
    >>> G5 = MixedMultiGraph(
    ...     undirected_edges=[(1, 2)],
    ...     directed_edges=[(2, 3), (3, 4)]
    ... )
    >>> G6 = MixedMultiGraph(
    ...     undirected_edges=[(5, 6)],
    ...     directed_edges=[(6, 7), (7, 8)]
    ... )
    >>> is_isomorphic(G5, G6)
    True
    
    Notes
    -----
    - Graph attributes are checked first for efficiency (early exit if they don't match).
    - NetworkX's efficient categorical matching functions are used internally.
    - Undirected edges are converted to bidirectional directed edges for matching.
    """
    # Check graph attributes first (early exit if they don't match)
    if graph_attrs:
        for attr in graph_attrs:
            if G1._directed.graph.get(attr) != G2._directed.graph.get(attr):
                return False
    
    # Convert to directed graphs
    nx_G1 = _to_digraph_for_isomorphism(G1)
    nx_G2 = _to_digraph_for_isomorphism(G2)
    
    # Create node match function
    if node_attrs:
        # Use NetworkX's efficient categorical_node_match
        # Default values are None, meaning nodes without the attribute match
        # only with nodes that also don't have it
        node_match = nx.isomorphism.categorical_node_match(
            node_attrs, [None] * len(node_attrs)
        )
    else:
        node_match = None
    
    # Create edge match function
    if edge_attrs:
        # For multi-graphs, use categorical_multiedge_match
        # Default values are None, meaning edges without the attribute match
        # only with edges that also don't have it
        edge_match = nx.isomorphism.categorical_multiedge_match(
            edge_attrs, [None] * len(edge_attrs)
        )
    else:
        edge_match = None
    
    # Use MultiDiGraphMatcher
    matcher = nx.isomorphism.MultiDiGraphMatcher(
        nx_G1, nx_G2,
        node_match=node_match,
        edge_match=edge_match
    )
    
    return matcher.is_isomorphic()


def _get_graph_invariant(graph: MixedMultiGraph) -> tuple[int, int, tuple[int, ...], tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    """
    Compute graph invariants for fast isomorphism candidate filtering.
    
    Returns a tuple of (num_nodes, num_edges, sorted_in_degrees, sorted_out_degrees, 
    sorted_undirected_degrees, sorted_edge_multiplicities). Isomorphic graphs must have 
    the same invariants (but not vice versa).
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to compute invariants for.
    
    Returns
    -------
    tuple[int, int, tuple[int, ...], tuple[int, ...], tuple[int, ...], tuple[int, ...]]
        Tuple of (num_nodes, num_edges, sorted_in_degrees, sorted_out_degrees, 
        sorted_undirected_degrees, sorted_edge_multiplicities). Edge multiplicities are 
        the number of parallel edges for each (u, v) pair, sorted.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> G = MixedMultiGraph(
    ...     directed_edges=[(1, 2), (1, 2)],
    ...     undirected_edges=[(2, 3)]
    ... )
    >>> inv = _get_graph_invariant(G)
    >>> inv[0]  # num_nodes
    3
    >>> inv[1]  # num_edges
    3
    """
    nodes = list(graph.nodes())
    num_nodes = len(nodes)
    num_edges = graph.number_of_edges()
    in_degrees = tuple(sorted(graph.indegree(v) for v in nodes))
    out_degrees = tuple(sorted(graph.outdegree(v) for v in nodes))
    undirected_degrees = tuple(sorted(graph._undirected.degree(v) for v in nodes))
    
    # Count edge multiplicities (number of parallel edges for each (u, v) pair)
    edge_multiplicities: list[int] = []
    seen_pairs: set[tuple[Any, Any]] = set()
    
    # Count directed edge multiplicities
    for u, v, _ in graph.directed_edges_iter(keys=True):
        if (u, v) not in seen_pairs:
            multiplicity = graph._directed.number_of_edges(u, v)
            edge_multiplicities.append(multiplicity)
            seen_pairs.add((u, v))
    
    # Count undirected edge multiplicities
    for u, v, _ in graph.undirected_edges_iter(keys=True):
        # Normalize edge for undirected (order doesn't matter)
        edge_key = MixedMultiGraph.normalize_undirected_edge(u, v)
        if edge_key not in seen_pairs:
            multiplicity = graph._undirected.number_of_edges(u, v)
            edge_multiplicities.append(multiplicity)
            seen_pairs.add(edge_key)
    
    sorted_multiplicities = tuple(sorted(edge_multiplicities))
    
    return (num_nodes, num_edges, in_degrees, out_degrees, undirected_degrees, sorted_multiplicities)

