"""
Network conversion module.

This module provides functions for converting between different graph representations
and semi-directed phylogenetic networks.
"""

from __future__ import annotations

from typing import Any, Literal, TypeVar

import networkx as nx

from .base import MixedPhyNetwork
from .sd_phynetwork import SemiDirectedPhyNetwork
from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.conversions import (
    graph_to_mixedmultigraph,
    multigraph_to_mixedmultigraph,
)
from ....utils.exceptions import PhyloZooValueError, PhyloZooTypeError
T = TypeVar('T')


def _sdnetwork_from_mmgraph(
    graph: MixedMultiGraph[T],
    network_type: Literal['semi-directed', 'mixed'] = 'semi-directed'
) -> SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]:
    """
    Internal helper to create a SemiDirectedPhyNetwork or MixedPhyNetwork from a MixedMultiGraph.

    Parameters
    ----------
    graph : MixedMultiGraph[T]
        The mixed multigraph to convert.
    network_type : Literal['semi-directed', 'mixed'], default='semi-directed'
        Type of network to create. 'semi-directed' creates a SemiDirectedPhyNetwork,
        'mixed' creates a MixedPhyNetwork.

    Returns
    -------
    SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]
        A new phylogenetic network with edges and labels from the graph.
    """
    # Extract directed edges
    directed_edges: list[dict[str, Any]] = []
    for u, v, key, data in graph.directed_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        if data:
            edge_dict.update(data)
        directed_edges.append(edge_dict)

    # Extract undirected edges
    undirected_edges: list[dict[str, Any]] = []
    for u, v, key, data in graph.undirected_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        if data:
            edge_dict.update(data)
        undirected_edges.append(edge_dict)

    # Extract nodes with their attributes
    nodes: list[tuple[T, dict[str, Any]]] = []
    for node in graph.nodes():
        attrs = graph._undirected.nodes[node]
        nodes.append((node, attrs.copy()))

    # Extract graph attributes
    graph_attributes = graph._directed.graph.copy()

    # Create and return new network
    if network_type == 'semi-directed':
        return SemiDirectedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            nodes=nodes if nodes else None,
            attributes=graph_attributes if graph_attributes else None
        )
    else:  # network_type == 'mixed'
        return MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            nodes=nodes if nodes else None,
            attributes=graph_attributes if graph_attributes else None
        )


def sdnetwork_from_graph(
    graph: nx.Graph | nx.MultiGraph | MixedMultiGraph[T],
    network_type: Literal['semi-directed', 'mixed'] = 'semi-directed'
) -> SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]:
    """
    Create a SemiDirectedPhyNetwork or MixedPhyNetwork from a NetworkX Graph, MultiGraph, or phylozoo
    MixedMultiGraph.
    
    For NetworkX graphs, all edges are treated as undirected edges. Edge attributes,
    node attributes, and graph-level attributes are preserved and passed through to
    the resulting network.
    
    Parameters
    ----------
    graph : nx.Graph | nx.MultiGraph | MixedMultiGraph[T]
        The graph to convert. Can be a NetworkX Graph, MultiGraph, or a
        MixedMultiGraph from the primitives module.
    network_type : Literal['semi-directed', 'mixed'], default='semi-directed'
        Type of network to create. 'semi-directed' creates a SemiDirectedPhyNetwork,
        'mixed' creates a MixedPhyNetwork.
    
    Returns
    -------
    SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]
        A new phylogenetic network with edges and labels from the graph.
    
    Raises
    ------
    PhyloZooValueError
        If the resulting network is invalid according to SemiDirectedPhyNetwork or
        MixedPhyNetwork validation rules (e.g., invalid node degrees, undirected
        cycles in semi-directed networks, etc.).
    
    Notes
    -----

    - **Edge attributes**: All edge attributes (e.g., `branch_length`, `bootstrap`,
      `gamma` for hybrid edges) are preserved and passed through to the network.
    - **Node attributes**: All node attributes are preserved. The `label` attribute
      is used for taxon labels on leaf nodes.
    - **Graph attributes**: Graph-level attributes are preserved and stored in the
      network's attributes dictionary.
    - **Validation**: The network is validated upon creation. If the graph structure
      does not meet network requirements (e.g., leaves must have no outgoing edges,
      internal nodes must have appropriate degrees, etc.), a ValueError is raised.
    
    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edge(0, 1, branch_length=0.5)
    >>> G.add_edge(0, 2, branch_length=0.3)
    >>> G.nodes[1]['label'] = 'A'
    >>> G.nodes[2]['label'] = 'B'
    >>> G.graph['source'] = 'test'
    >>> net = sdnetwork_from_graph(G)
    >>> isinstance(net, SemiDirectedPhyNetwork)
    True
    >>> net.get_label(1)
    'A'
    >>> net.get_branch_length(0, 1)
    0.5
    >>> net.get_network_attribute('source')
    'test'
    """
    # Convert NetworkX graph to MixedMultiGraph if needed
    if isinstance(graph, (nx.Graph, nx.MultiGraph)):
        if isinstance(graph, nx.MultiGraph):
            mmgraph = multigraph_to_mixedmultigraph(graph)
        else:
            mmgraph = graph_to_mixedmultigraph(graph)
    elif isinstance(graph, MixedMultiGraph):
        mmgraph = graph
    else:
        raise PhyloZooTypeError(
            f"Expected nx.Graph, nx.MultiGraph, or MixedMultiGraph, "
            f"got {type(graph)}"
        )
    
    # Convert MixedMultiGraph to network
    return _sdnetwork_from_mmgraph(mmgraph, network_type=network_type)

