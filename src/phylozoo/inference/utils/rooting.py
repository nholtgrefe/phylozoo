"""
Network rooting utilities.

This module provides functions for rooting phylogenetic networks.
"""

from __future__ import annotations

from typing import Any

from ...core.network.dnetwork import DirectedPhyNetwork
from ...core.network.sdnetwork import SemiDirectedPhyNetwork
from ...core.network.sdnetwork.derivations import to_d_network
from ...utils.exceptions import PhyloZooValueError

def root_at_outgroup(
    network: SemiDirectedPhyNetwork,
    outgroup: str,
) -> DirectedPhyNetwork:
    """
    Root a semi-directed network at the edge leading to the specified outgroup taxon.
    
    This function finds the edge incident to the leaf node with the given taxon label,
    then roots the network at that edge using `to_d_network`.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to root.
    outgroup : str
        The taxon label of the outgroup leaf.
    
    Returns
    -------
    DirectedPhyNetwork
        A directed phylogenetic network rooted at the edge leading to the outgroup.
    
    Raises
    ------
    PhyloZooValueError
        If the outgroup taxon is not found in the network.
        If no edge is found incident to the outgroup leaf.
        If the edge is not a valid root location.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.inference.utils import root_at_outgroup
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> d_net = root_at_outgroup(net, 'A')
    >>> d_net.root_node is not None
    True
    """
    # Get the node ID for the outgroup taxon
    outgroup_node = network.get_node_id(outgroup)
    if outgroup_node is None:
        raise PhyloZooValueError(f"Outgroup taxon '{outgroup}' not found in network")
    
    # Find an edge incident to the outgroup leaf node
    # Check undirected edges first (most common case)
    incident_edges: list[tuple[Any, Any, int]] = []
    
    for u, v, key in network._graph.incident_undirected_edges(outgroup_node, keys=True):
        # The edge is (u, v, key) where one of u or v is outgroup_node
        # We want the edge as (parent, outgroup_node, key) for rooting
        if u == outgroup_node:
            incident_edges.append((v, u, key))
        else:
            incident_edges.append((u, v, key))
    
    if not incident_edges:
        raise PhyloZooValueError(
            f"No edge found incident to outgroup leaf node {outgroup_node} "
            f"(taxon '{outgroup}')"
        )
    
    # Use the first incident edge as the root location
    # For undirected edges, we can root on the edge
    root_location = incident_edges[0]
    
    # Convert to directed network
    return to_d_network(network, root_location=root_location)
