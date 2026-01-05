"""
Isomorphism module for mixed multi-graphs.

This module provides functions for checking graph isomorphism between
MixedMultiGraph instances.
"""

from __future__ import annotations

from typing import TypeVar

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

