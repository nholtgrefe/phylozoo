"""
Conversion functions for MixedMultiGraph.

This module provides functions for converting NetworkX graphs to MixedMultiGraph.
"""

from typing import TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from ..d_multigraph import DirectedMultiGraph
    from . import MixedMultiGraph
else:
    from . import MixedMultiGraph


def graph_to_mixedmultigraph(graph: nx.Graph) -> 'MixedMultiGraph':
    """
    Create a MixedMultiGraph from a NetworkX Graph.
    
    All edges from the Graph are added as undirected edges.
    
    Parameters
    ----------
    graph : nx.Graph
        NetworkX Graph to convert.
    
    Returns
    -------
    MixedMultiGraph
        New MixedMultiGraph instance with all edges as undirected.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edge(1, 2, weight=5.0)
    >>> G.add_edge(2, 3)
    >>> M = graph_to_mixedmultigraph(G)
    >>> M.number_of_edges()
    2
    """
    mg = MixedMultiGraph()
    # Add all nodes with attributes
    for node, data in graph.nodes(data=True):
        mg.add_node(node, **data)
    # Add all edges with attributes
    for u, v, data in graph.edges(data=True):
        mg.add_undirected_edge(u, v, **data)
    return mg


def multigraph_to_mixedmultigraph(graph: nx.MultiGraph) -> 'MixedMultiGraph':
    """
    Create a MixedMultiGraph from a NetworkX MultiGraph.
    
    All edges from the MultiGraph are added as undirected edges, preserving
    parallel edges and their keys.
    
    Parameters
    ----------
    graph : nx.MultiGraph
        NetworkX MultiGraph to convert.
    
    Returns
    -------
    MixedMultiGraph
        New MixedMultiGraph instance with all edges as undirected.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.MultiGraph()
    >>> G.add_edge(1, 2, key=0, weight=1.0)
    0
    >>> G.add_edge(1, 2, key=1, weight=2.0)
    1
    >>> M = multigraph_to_mixedmultigraph(G)
    >>> M.number_of_edges()
    2
    """
    mg = MixedMultiGraph()
    # Add all nodes with attributes
    for node, data in graph.nodes(data=True):
        mg.add_node(node, **data)
    # Add all edges with keys and attributes
    for u, v, key, data in graph.edges(keys=True, data=True):
        mg.add_undirected_edge(u, v, key=key, **data)
    return mg


def multidigraph_to_mixedmultigraph(graph: nx.MultiDiGraph) -> 'MixedMultiGraph':
    """
    Create a MixedMultiGraph from a NetworkX MultiDiGraph.
    
    All edges from the MultiDiGraph are added as directed edges, preserving
    parallel edges and their keys.
    
    Parameters
    ----------
    graph : nx.MultiDiGraph
        NetworkX MultiDiGraph to convert.
    
    Returns
    -------
    MixedMultiGraph
        New MixedMultiGraph instance with all edges as directed.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.MultiDiGraph()
    >>> G.add_edge(1, 2, key=0, weight=1.0)
    0
    >>> G.add_edge(1, 2, key=1, weight=2.0)
    1
    >>> M = multidigraph_to_mixedmultigraph(G)
    >>> M._directed.number_of_edges()
    2
    """
    mg = MixedMultiGraph()
    # Add all nodes with attributes
    for node, data in graph.nodes(data=True):
        mg.add_node(node, **data)
    # Add all edges directly to preserve keys (bypass mutual exclusivity checks)
    for u, v, key, data in graph.edges(keys=True, data=True):
        mg._directed.add_edge(u, v, key=key, **data)
        mg._combined.add_edge(u, v, key=key, **data)
    return mg


def directedmultigraph_to_mixedmultigraph(graph: 'DirectedMultiGraph') -> 'MixedMultiGraph':
    """
    Create a MixedMultiGraph from a DirectedMultiGraph.
    
    All edges from the DirectedMultiGraph are added as directed edges,
    preserving parallel edges and their keys.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        DirectedMultiGraph instance to convert.
    
    Returns
    -------
    MixedMultiGraph
        New MixedMultiGraph instance with all edges as directed.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2, weight=1.0)
    0
    >>> M = directedmultigraph_to_mixedmultigraph(G)
    >>> M._directed.number_of_edges()
    1
    """
    mg = MixedMultiGraph()
    # Add all nodes with attributes
    for node, data in graph.nodes(data=True):
        mg.add_node(node, **data)
    # Add all edges directly to preserve keys (bypass mutual exclusivity checks)
    for u, v, key, data in graph.edges(keys=True, data=True):
        mg._directed.add_edge(u, v, key=key, **data)
        mg._combined.add_edge(u, v, key=key, **data)
    return mg

