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
from ....utils.exceptions.utils import warn_on_keyword, warn_on_none_value
from ...primitives.d_multigraph import DirectedMultiGraph
from ...primitives.d_multigraph.conversions import (
    digraph_to_directedmultigraph,
    multidigraph_to_directedmultigraph,
)

T = TypeVar("T")


def _dnetwork_from_dmgraph(
    graph: DirectedMultiGraph[T], *, copy: bool = True
) -> DirectedPhyNetwork[T]:
    """
    Internal helper to create a DirectedPhyNetwork from a DirectedMultiGraph.

    Parameters
    ----------
    graph : DirectedMultiGraph[T]
        The directed multigraph to convert.
    copy : bool, optional
        If True, the network is rebuilt from the graph's edges and nodes and never
        shares state with ``graph``. If False, the network adopts ``graph`` itself:
        the same label checks, leaf auto-labelling and :meth:`validate` run, but the
        graph is not constructed a second time. Only pass False for a graph the
        caller owns and will not touch again. By default True.

    Returns
    -------
    DirectedPhyNetwork[T]
        A new directed phylogenetic network with edges and labels from the graph.
    """
    if not copy:
        return _adopt_dmgraph(graph)

    # Extract edges
    edges: list[dict[str, Any]] = []
    for u, v, key, data in graph.edges(keys=True, data=True):
        edge_dict: dict[str, Any] = {"u": u, "v": v}
        if key != 0:
            edge_dict["key"] = key
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
        edges=edges, nodes=nodes, attributes=graph_attributes if graph_attributes else None
    )


def _adopt_dmgraph(graph: DirectedMultiGraph[T]) -> DirectedPhyNetwork[T]:
    """
    Build a DirectedPhyNetwork around an existing DirectedMultiGraph without rebuilding it.

    Performs the constructor's steps after graph construction, on the graph as given:
    the per-node checks ``add_node`` would make, label registration with the same
    string/uniqueness validation, leaf auto-labelling, and :meth:`validate`.

    Parameters
    ----------
    graph : DirectedMultiGraph[T]
        The graph the network takes ownership of.

    Returns
    -------
    DirectedPhyNetwork[T]
        The network wrapping ``graph``.
    """
    network = DirectedPhyNetwork.__new__(DirectedPhyNetwork)
    network._graph = graph
    network._node_to_label = {}
    network._label_to_node = {}
    for node, attrs in graph._graph.nodes(data=True):
        warn_on_keyword(node, "Node id")
        for attr_name, attr_value in attrs.items():
            warn_on_keyword(attr_name, "Attribute name")
            warn_on_none_value(attr_value, f"Attribute '{attr_name}'")
        if "label" in attrs:
            network._add_label_to_dicts(node, attrs["label"])
    network._auto_label_unlabeled_leaves()
    network.validate()
    return network


def dnetwork_from_graph(
    graph: nx.DiGraph | nx.MultiDiGraph | DirectedMultiGraph[T],
    copy: bool = True,
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
    copy : bool, optional
        Only relevant when ``graph`` is a DirectedMultiGraph. If True, the network
        is built from a fresh copy of the graph's contents. If False, the network
        takes ownership of ``graph`` itself and the caller must not use or modify
        it afterwards; this skips one full graph construction. Validation is the
        same either way. NetworkX inputs are always converted into a new graph.
        By default True.

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
        # The converted graph is brand new and owned here, so it can be adopted.
        copy = False
    elif isinstance(graph, DirectedMultiGraph):
        dmgraph = graph
    else:
        raise PhyloZooTypeError(
            f"Expected nx.DiGraph, nx.MultiDiGraph, or DirectedMultiGraph, " f"got {type(graph)}"
        )

    # Convert DirectedMultiGraph to network
    return _dnetwork_from_dmgraph(dmgraph, copy=copy)
