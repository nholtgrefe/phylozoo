"""
Isomorphism module for semi-directed and mixed phylogenetic networks.

This module provides functions for checking network isomorphism between
SemiDirectedPhyNetwork and MixedPhyNetwork instances.
"""

from __future__ import annotations

from typing import TypeVar

from ...primitives.m_multigraph.isomorphism import is_isomorphic as mm_is_isomorphic
from .base import MixedPhyNetwork

T = TypeVar('T')


def is_isomorphic(
    net1: MixedPhyNetwork[T],
    net2: MixedPhyNetwork[T],
    node_attrs: list[str] | None = None,
    edge_attrs: list[str] | None = None,
    graph_attrs: list[str] | None = None,
) -> bool:
    """
    Check if two semi-directed or mixed phylogenetic networks are isomorphic.
    
    Two networks are isomorphic if there exists a bijection between their node sets
    that preserves adjacency, edge direction, parallel edges, and node labels.
    Labels are always checked (non-optional), and additional node, edge, and graph
    attributes can be specified.
    
    Parameters
    ----------
    net1 : MixedPhyNetwork
        First phylogenetic network (SemiDirectedPhyNetwork or MixedPhyNetwork).
    net2 : MixedPhyNetwork
        Second phylogenetic network (SemiDirectedPhyNetwork or MixedPhyNetwork).
    node_attrs : list[str] | None, optional
        List of additional node attribute names to match (beyond 'label').
        If None, only 'label' is checked. By default None.
    edge_attrs : list[str] | None, optional
        List of edge attribute names to match. If None, edge attributes are ignored.
        By default None.
    graph_attrs : list[str] | None, optional
        List of graph-level attribute names to match. If None, graph attributes are
        ignored. By default None.
    
    Returns
    -------
    bool
        True if the networks are isomorphic, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net1 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> net2 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(4, 5), (4, 6)],
    ...     nodes=[(5, {'label': 'A'}), (6, {'label': 'B'})]
    ... )
    >>> is_isomorphic(net1, net2)
    True
    >>> # Different labels: not isomorphic
    >>> net3 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(4, 5), (4, 6)],
    ...     nodes=[(5, {'label': 'A'}), (6, {'label': 'C'})]
    ... )
    >>> is_isomorphic(net1, net3)
    False
    >>> # With additional node attributes
    >>> net4 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A', 'type': 'leaf'}), (2, {'label': 'B', 'type': 'leaf'})]
    ... )
    >>> net5 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(4, 5), (4, 6)],
    ...     nodes=[(5, {'label': 'A', 'type': 'leaf'}), (6, {'label': 'B', 'type': 'leaf'})]
    ... )
    >>> is_isomorphic(net4, net5, node_attrs=['type'])
    True
    >>> # With edge attributes
    >>> net6 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
    ...     nodes=[(1, {'label': 'A'})]
    ... )
    >>> net7 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[{'u': 4, 'v': 5, 'branch_length': 0.5}],
    ...     nodes=[(5, {'label': 'A'})]
    ... )
    >>> is_isomorphic(net6, net7, edge_attrs=['branch_length'])
    True
    
    Notes
    -----

    - Labels are always checked (non-optional) to ensure networks with different
      taxon labels are not considered isomorphic.
    - The function uses the underlying MixedMultiGraph isomorphism checking,
      which converts undirected edges to bidirectional directed edges for matching.
    - For node attributes, if a node doesn't have an attribute, it only matches with
      nodes that also don't have that attribute (None matches None).
    """
    # Always include 'label' in node attributes
    if node_attrs is None:
        node_attrs_list = ['label']
    else:
        # Ensure 'label' is included (avoid duplicates)
        node_attrs_list = ['label'] + [attr for attr in node_attrs if attr != 'label']
    
    # Call the underlying graph isomorphism function
    return mm_is_isomorphic(
        net1._graph,
        net2._graph,
        node_attrs=node_attrs_list,
        edge_attrs=edge_attrs,
        graph_attrs=graph_attrs,
    )

