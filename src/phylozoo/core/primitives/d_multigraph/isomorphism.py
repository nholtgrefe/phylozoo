"""
Isomorphism module for directed multi-graphs.

This module provides functions for checking graph isomorphism between
DirectedMultiGraph instances.
"""

from __future__ import annotations

from typing import TypeVar

import networkx as nx

from .base import DirectedMultiGraph

T = TypeVar('T')


def is_isomorphic(
    G1: DirectedMultiGraph[T],
    G2: DirectedMultiGraph[T],
    node_attrs: list[str] | None = None,
    edge_attrs: list[str] | None = None,
    graph_attrs: list[str] | None = None,
) -> bool:
    """
    Check if two directed multi-graphs are isomorphic.
    
    Two graphs are isomorphic if there exists a bijection between their node sets
    that preserves adjacency, edge direction, and parallel edges. Optionally,
    node attributes, edge attributes, and graph-level attributes can be required
    to match as well.
    
    Parameters
    ----------
    G1 : DirectedMultiGraph
        First directed multi-graph.
    G2 : DirectedMultiGraph
        Second directed multi-graph.
    node_attrs : list[str] | None, optional
        List of node attribute names to match. If None, node attributes are ignored.
        Nodes must have matching values for all specified attributes to be considered
        isomorphic. If a node doesn't have an attribute, it matches only with nodes
        that also don't have that attribute. By default None.
    edge_attrs : list[str] | None, optional
        List of edge attribute names to match. If None, edge attributes are ignored.
        Edges must have matching values for all specified attributes. If an edge
        doesn't have an attribute, it matches only with nodes that also don't have 
        that attribute. By default None.
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
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
    >>> G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
    >>> is_isomorphic(G1, G2)
    True
    >>> # With node labels: graphs are isomorphic if labels match
    >>> G1.add_node(1, label='A')
    >>> G1.add_node(2, label='B')
    >>> G2.add_node(4, label='A')
    >>> G2.add_node(5, label='B')
    >>> is_isomorphic(G1, G2, node_attrs=['label'])
    True
    >>> # Different labels: not isomorphic
    >>> G2.add_node(4, label='C')  # Update label
    >>> is_isomorphic(G1, G2, node_attrs=['label'])
    False
    >>> # Multiple node attributes
    >>> G1.add_node(1, type='leaf')
    >>> G2.add_node(4, type='leaf')
    >>> is_isomorphic(G1, G2, node_attrs=['label', 'type'])
    False  # Still false because label 'C' != 'A'
    >>> # Edge attributes
    >>> G3 = DirectedMultiGraph(edges=[{'u': 1, 'v': 2, 'weight': 1.0}])
    >>> G4 = DirectedMultiGraph(edges=[{'u': 3, 'v': 4, 'weight': 1.0}])
    >>> is_isomorphic(G3, G4, edge_attrs=['weight'])
    True
    >>> # Graph attributes
    >>> G5 = DirectedMultiGraph(attributes={'version': '1.0'})
    >>> G6 = DirectedMultiGraph(attributes={'version': '1.0'})
    >>> is_isomorphic(G5, G6, graph_attrs=['version'])
    True
    >>> # Parallel edges are preserved
    >>> G7 = DirectedMultiGraph(edges=[(1, 2), (1, 2)])
    >>> G8 = DirectedMultiGraph(edges=[(3, 4), (3, 4)])
    >>> is_isomorphic(G7, G8)
    True
    >>> # Different number of parallel edges: not isomorphic
    >>> G9 = DirectedMultiGraph(edges=[(1, 2)])
    >>> is_isomorphic(G7, G9)
    False
    
    Notes
    -----
    - Graph attributes are checked first for efficiency (early exit if they don't match).
    - NetworkX's efficient categorical matching functions are used internally.
    """
    # Get the underlying NetworkX graphs
    nx_G1 = G1._graph
    nx_G2 = G2._graph
    
    # Check graph attributes first (early exit if they don't match)
    if graph_attrs:
        for attr in graph_attrs:
            if nx_G1.graph.get(attr) != nx_G2.graph.get(attr):
                return False
    
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
    
    # Use MultiDiGraphMatcher with matching functions
    matcher = nx.isomorphism.MultiDiGraphMatcher(
        nx_G1, nx_G2,
        node_match=node_match,
        edge_match=edge_match
    )
    
    return matcher.is_isomorphic()

