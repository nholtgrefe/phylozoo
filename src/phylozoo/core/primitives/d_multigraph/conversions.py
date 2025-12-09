"""
Conversion functions for DirectedMultiGraph.

This module provides functions for converting NetworkX graphs to DirectedMultiGraph.
"""

import networkx as nx

from . import DirectedMultiGraph


def digraph_to_directedmultigraph(graph: nx.DiGraph) -> 'DirectedMultiGraph':
    """
    Create a DirectedMultiGraph from a NetworkX DiGraph.
    
    All edges from the DiGraph are added as directed edges.
    
    Parameters
    ----------
    graph : nx.DiGraph
        NetworkX DiGraph to convert.
    
    Returns
    -------
    DirectedMultiGraph
        New DirectedMultiGraph instance with all edges as directed.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.DiGraph()
    >>> G.add_edge(1, 2, weight=5.0)
    >>> G.add_edge(2, 3)
    >>> M = digraph_to_directedmultigraph(G)
    >>> M.number_of_edges()
    2
    """
    dmg = DirectedMultiGraph()
    # Add all nodes with attributes
    for node, data in graph.nodes(data=True):
        dmg.add_node(node, **data)
    # Add all edges with attributes
    for u, v, data in graph.edges(data=True):
        dmg.add_edge(u, v, **data)
    return dmg


def multidigraph_to_directedmultigraph(graph: nx.MultiDiGraph) -> 'DirectedMultiGraph':
    """
    Create a DirectedMultiGraph from a NetworkX MultiDiGraph.
    
    All edges from the MultiDiGraph are added as directed edges, preserving
    parallel edges and their keys.
    
    Parameters
    ----------
    graph : nx.MultiDiGraph
        NetworkX MultiDiGraph to convert.
    
    Returns
    -------
    DirectedMultiGraph
        New DirectedMultiGraph instance with all edges as directed.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.MultiDiGraph()
    >>> G.add_edge(1, 2, key=0, weight=1.0)
    0
    >>> G.add_edge(1, 2, key=1, weight=2.0)
    1
    >>> M = multidigraph_to_directedmultigraph(G)
    >>> M.number_of_edges()
    2
    """
    dmg = DirectedMultiGraph()
    # Add all nodes with attributes
    for node, data in graph.nodes(data=True):
        dmg.add_node(node, **data)
    # Add all edges with keys and attributes
    for u, v, key, data in graph.edges(keys=True, data=True):
        dmg.add_edge(u, v, key=key, **data)
    return dmg

