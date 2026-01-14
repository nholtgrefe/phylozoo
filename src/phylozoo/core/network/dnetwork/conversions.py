"""
Network conversion module.

This module provides functions for converting between different graph representations
and directed phylogenetic networks.
"""

from __future__ import annotations

from typing import Any, TypeVar

import networkx as nx

from ....utils.exceptions import PhyloZooTypeError
from .base import DirectedPhyNetwork
from ...primitives.d_multigraph import DirectedMultiGraph
from ...primitives.d_multigraph.conversions import (
    digraph_to_directedmultigraph,
    multidigraph_to_directedmultigraph,
)

T = TypeVar('T')


def _dnetwork_from_dmgraph(graph: DirectedMultiGraph[T]) -> DirectedPhyNetwork[T]:
    """
    Internal helper to create a DirectedPhyNetwork from a DirectedMultiGraph.

    Parameters
    ----------
    graph : DirectedMultiGraph[T]
        The directed multigraph to convert.

    Returns
    -------
    DirectedPhyNetwork[T]
        A new directed phylogenetic network with edges and labels from the graph.
    """
    # Extract edges
    edges: list[dict[str, Any]] = []
    for u, v, key, data in graph.edges(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        if data:
            edge_dict.update(data)
        edges.append(edge_dict)

    # Extract nodes with their attributes
    nodes: list[tuple[T, dict[str, Any]]] = []
    for node in graph.nodes():
        attrs = graph._graph.nodes[node]
        nodes.append((node, attrs.copy()))

    # Extract graph attributes
    graph_attributes = graph._graph.graph.copy()

    # Create and return new network
    return DirectedPhyNetwork(
        edges=edges,
        nodes=nodes,
        attributes=graph_attributes if graph_attributes else None
    )


def dnetwork_from_graph(
    graph: nx.DiGraph | nx.MultiDiGraph | DirectedMultiGraph[T],
) -> DirectedPhyNetwork[T]:
    """
    Create a DirectedPhyNetwork from a NetworkX DiGraph, MultiDiGraph, or phylozoo
    DirectedMultiGraph.
    
    All edges from the input graph are treated as directed edges. Edge attributes,
    node attributes, and graph-level attributes are preserved and passed through to
    the resulting network.
    
    Parameters
    ----------
    graph : nx.DiGraph | nx.MultiDiGraph | DirectedMultiGraph[T]
        The graph to convert. Can be a NetworkX DiGraph, MultiDiGraph, or a
        DirectedMultiGraph from the primitives module.
    
    Returns
    -------
    DirectedPhyNetwork[T]
        A new directed phylogenetic network with edges and labels from the graph.
    
    Raises
    ------
    PhyloZooTypeError
        If graph is not one of the supported types (nx.DiGraph, nx.MultiDiGraph, or DirectedMultiGraph).
    PhyloZooValueError
        If the resulting network is invalid according to DirectedPhyNetwork validation
        rules (e.g., not a DAG, invalid node degrees, etc.).
    
    Notes
    -----
    - **Edge attributes**: All edge attributes (e.g., `branch_length`, `gamma`,
      `bootstrap`) are preserved and passed through to the network.
    - **Node attributes**: All node attributes are preserved. The `label` attribute
      is used for taxon labels on leaf nodes.
    - **Graph attributes**: Graph-level attributes are preserved and stored in the
      network's attributes dictionary.
    - **Validation**: The network is validated upon creation. If the graph structure
      does not meet DirectedPhyNetwork requirements (e.g., must be a DAG, leaves
      must have in-degree 1, etc.), a ValueError is raised.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.DiGraph()
    >>> G.add_edge(0, 1, branch_length=0.5)
    >>> G.add_edge(0, 2, branch_length=0.3)
    >>> G.nodes[1]['label'] = 'A'
    >>> G.nodes[2]['label'] = 'B'
    >>> G.graph['source'] = 'test'
    >>> net = dnetwork_from_graph(G)
    >>> isinstance(net, DirectedPhyNetwork)
    True
    >>> net.get_label(1)
    'A'
    >>> net.get_branch_length(0, 1)
    0.5
    >>> net.get_network_attribute('source')
    'test'
    """
    # Convert NetworkX graph to DirectedMultiGraph if needed
    if isinstance(graph, (nx.DiGraph, nx.MultiDiGraph)):
        if isinstance(graph, nx.MultiDiGraph):
            dmgraph = multidigraph_to_directedmultigraph(graph)
        else:
            dmgraph = digraph_to_directedmultigraph(graph)
    elif isinstance(graph, DirectedMultiGraph):
        dmgraph = graph
    else:
        raise PhyloZooTypeError(
            f"Expected nx.DiGraph, nx.MultiDiGraph, or DirectedMultiGraph, "
            f"got {type(graph)}"
        )
    
    # Convert DirectedMultiGraph to network
    return _dnetwork_from_dmgraph(dmgraph)

