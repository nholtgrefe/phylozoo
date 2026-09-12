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
from ....utils.exceptions.utils import warn_on_keyword, warn_on_none_value
from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.conversions import (
    graph_to_mixedmultigraph,
    multigraph_to_mixedmultigraph,
)
from ....utils.exceptions import PhyloZooTypeError

T = TypeVar("T")


def _sdnetwork_from_mmgraph(
    graph: MixedMultiGraph[T],
    network_type: Literal["semi-directed", "mixed"] = "semi-directed",
    *,
    copy: bool = True,
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

    copy : bool, optional
        If True, the network is rebuilt from the graph's edges and nodes and never
        shares state with ``graph``. If False, the network adopts ``graph`` itself:
        the same label checks, leaf auto-labelling and :meth:`validate` run, but the
        graph is not constructed a second time. Only pass False for a graph the
        caller owns and will not touch again. By default True.

    Returns
    -------
    SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]
        A new phylogenetic network with edges and labels from the graph.
    """
    if not copy:
        return _adopt_mmgraph(graph, network_type)

    # Extract directed edges
    directed_edges: list[dict[str, Any]] = []
    for u, v, key, data in graph.directed_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {"u": u, "v": v}
        if key != 0:
            edge_dict["key"] = key
        if data:
            edge_dict.update(data)
        directed_edges.append(edge_dict)

    # Extract undirected edges
    undirected_edges: list[dict[str, Any]] = []
    for u, v, key, data in graph.undirected_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {"u": u, "v": v}
        if key != 0:
            edge_dict["key"] = key
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
    if network_type == "semi-directed":
        return SemiDirectedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            nodes=nodes if nodes else None,
            attributes=graph_attributes if graph_attributes else None,
        )
    else:  # network_type == 'mixed'
        return MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            nodes=nodes if nodes else None,
            attributes=graph_attributes if graph_attributes else None,
        )


def _adopt_mmgraph(
    graph: MixedMultiGraph[T],
    network_type: Literal["semi-directed", "mixed"],
) -> SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]:
    """
    Build a network around an existing MixedMultiGraph without rebuilding it.

    Performs the constructor's steps after graph construction, on the graph as
    given: the per-node checks ``add_node`` would make, the mirroring of each node's
    undirected-side attributes onto the directed side that re-adding the nodes used
    to do, label registration with the same string/uniqueness validation, leaf
    auto-labelling, and :meth:`validate`.

    Parameters
    ----------
    graph : MixedMultiGraph[T]
        The graph the network takes ownership of.
    network_type : Literal['semi-directed', 'mixed']
        Which network class to build.

    Returns
    -------
    SemiDirectedPhyNetwork[T] | MixedPhyNetwork[T]
        The network wrapping ``graph``.
    """
    cls = SemiDirectedPhyNetwork if network_type == "semi-directed" else MixedPhyNetwork
    network = cls.__new__(cls)
    network._graph = graph
    network._node_to_label = {}
    network._label_to_node = {}
    undirected_nodes = graph._undirected.nodes
    directed = graph._directed
    for node in graph.nodes():
        attrs = undirected_nodes[node]
        warn_on_keyword(node, "Node id")
        for attr_name, attr_value in attrs.items():
            warn_on_keyword(attr_name, "Attribute name")
            warn_on_none_value(attr_value, f"Attribute '{attr_name}'")
        # add_node(node, **attrs) set these on both sub-graphs; keep that invariant.
        directed.add_node(node, **attrs)
        if "label" in attrs:
            network._add_label_to_dicts(node, attrs["label"])
    graph._combined_cache = None
    network._auto_label_unlabeled_leaves()
    network.validate()
    return network


def sdnetwork_from_graph(
    graph: nx.Graph | nx.MultiGraph | MixedMultiGraph[T],
    network_type: Literal["semi-directed", "mixed"] = "semi-directed",
    copy: bool = True,
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
    copy : bool, optional
        Only relevant when ``graph`` is a MixedMultiGraph. If True, the network is
        built from a fresh copy of the graph's contents. If False, the network takes
        ownership of ``graph`` itself and the caller must not use or modify it
        afterwards; this skips one full graph construction. Validation is the same
        either way. NetworkX inputs are always converted into a new graph.
        By default True.

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
    >>> G.add_edge(0, 3, branch_length=0.2)
    >>> G.nodes[1]['label'] = 'A'
    >>> G.nodes[2]['label'] = 'B'
    >>> G.nodes[3]['label'] = 'C'
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
        # The converted graph is brand new and owned here, so it can be adopted.
        copy = False
    elif isinstance(graph, MixedMultiGraph):
        mmgraph = graph
    else:
        raise PhyloZooTypeError(
            f"Expected nx.Graph, nx.MultiGraph, or MixedMultiGraph, " f"got {type(graph)}"
        )

    # Convert MixedMultiGraph to network
    return _sdnetwork_from_mmgraph(mmgraph, network_type=network_type, copy=copy)
