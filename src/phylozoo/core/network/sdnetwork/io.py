"""
Network I/O module.

This module provides functions for reading and writing semi-directed and mixed phylogenetic networks
to/from various file formats (e.g., Newick, Nexus, eNewick, etc.).
"""

from __future__ import annotations

from typing import Any, Literal, TypeVar

from .base import MixedPhyNetwork
from .sd_phynetwork import SemiDirectedPhyNetwork
from ...primitives.m_multigraph import MixedMultiGraph

T = TypeVar('T')


def sdnetwork_from_mmgraph(
    graph: MixedMultiGraph[T],
    network_type: Literal['semi-directed', 'mixed'] = 'semi-directed'
) -> SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]:
    """
    Create a SemiDirectedPhyNetwork or MixedPhyNetwork from a MixedMultiGraph.

    Extracts directed edges, undirected edges, and node labels from the graph
    to create a new network. Node labels are extracted from the 'label' node attribute.

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

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> mmgraph = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
    >>> mmgraph.add_node(3, label='A')
    >>> net = sdnetwork_from_mmgraph(mmgraph, network_type='semi-directed')
    >>> isinstance(net, SemiDirectedPhyNetwork)
    True
    >>> net.get_label(3)
    'A'
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

