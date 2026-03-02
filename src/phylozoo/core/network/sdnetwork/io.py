"""
Network I/O module.

This module registers format handlers for SemiDirectedPhyNetwork with FormatRegistry
for use with the IOMixin system.

All format handlers are defined here:
- eNewick format: Uses to_d_network and to_sd_network for conversion
- PhyloZoo-DOT format: Uses MixedMultiGraph's phylozoo-dot I/O functions
"""

from __future__ import annotations

from typing import Any, TYPE_CHECKING

from phylozoo.utils.io import FormatRegistry

if TYPE_CHECKING:
    from .sd_phynetwork import SemiDirectedPhyNetwork


def to_newick(sd_network: 'SemiDirectedPhyNetwork', **kwargs: Any) -> str:
    """
    Convert a SemiDirectedPhyNetwork to an eNewick format string.
    
    This function converts the semi-directed network to a directed network
    using to_d_network, then converts the directed network to eNewick format.
    
    Parameters
    ----------
    sd_network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to convert.
    **kwargs
        Additional arguments to pass to to_d_network and to_enewick. Use
        root_location to specify RootLocation when converting to directed network
        (if None, a default location is chosen). Other arguments are passed to
        to_enewick.
    
    Returns
    -------
    str
        The eNewick format string representation of the network.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.io import to_newick
    >>> 
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> newick_str = to_newick(net)
    >>> ';' in newick_str
    True
    >>> 'A' in newick_str
    True
    >>> 'B' in newick_str
    True
    
    Notes
    -----
    The conversion process:
    1. Convert SemiDirectedPhyNetwork to DirectedPhyNetwork using to_d_network
    2. Convert DirectedPhyNetwork to eNewick string using to_enewick
    
    If no root_location is provided, a default location is chosen automatically.
    """
    # Local imports to avoid circular dependencies
    from .derivations import to_d_network
    from ..dnetwork._enewick import to_enewick
    
    # Extract root_location from kwargs if provided
    root_location = kwargs.pop('root_location', None)
    
    # Convert to directed network
    d_network = to_d_network(sd_network, root_location=root_location)
    
    # Convert to eNewick string
    return to_enewick(d_network, **kwargs)


def from_newick(newick_string: str, **kwargs: Any) -> 'SemiDirectedPhyNetwork':
    """
    Parse an eNewick format string and create a SemiDirectedPhyNetwork.
    
    This function parses the eNewick string into a DirectedPhyNetwork,
    then converts it to a SemiDirectedPhyNetwork using to_sd_network.
    
    Parameters
    ----------
    newick_string : str
        eNewick format string containing network data.
    **kwargs
        Additional arguments to pass to from_enewick (currently unused, for compatibility).
    
    Returns
    -------
    SemiDirectedPhyNetwork
        Parsed semi-directed phylogenetic network.
    
    Raises
    ------
    ENewickParseError
        If the eNewick string is malformed or cannot be parsed, or if the resulting
        network structure is invalid for SemiDirectedPhyNetwork.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork.io import from_newick
    >>> 
    >>> newick_str = "((A,B),C);"
    >>> net = from_newick(newick_str)
    >>> net.number_of_nodes()
    4
    >>> 'A' in net.taxa
    True
    
    Notes
    -----
    The conversion process:
    1. Parse eNewick string to DirectedPhyNetwork using from_enewick
    2. Convert DirectedPhyNetwork to SemiDirectedPhyNetwork using to_sd_network
    """
    # Local imports to avoid circular dependencies
    from ..dnetwork._enewick import from_enewick
    from ..dnetwork.derivations import to_sd_network
    
    # Parse eNewick string to directed network
    d_network = from_enewick(newick_string, **kwargs)
    
    # Convert to semi-directed network
    return to_sd_network(d_network)


def to_phylozoo_dot(sd_network: 'SemiDirectedPhyNetwork', **kwargs: Any) -> str:
    """
    Convert a SemiDirectedPhyNetwork to a PhyloZoo-DOT format string.
    
    This function delegates to MixedMultiGraph's to_phylozoo_dot function,
    using the underlying _graph attribute.
    
    Parameters
    ----------
    sd_network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to convert.
    **kwargs
        Additional arguments to pass to the underlying to_phylozoo_dot:
        - graph_name (str): Optional name for the graph (default: 'graph').
    
    Returns
    -------
    str
        The PhyloZoo-DOT format string representation of the network.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.io import to_phylozoo_dot
    >>> 
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> dot_str = to_phylozoo_dot(net)
    >>> 'graph' in dot_str
    True
    >>> '--' in dot_str
    True
    
    Notes
    -----
    The PhyloZoo-DOT format uses:
    - `--` for undirected edges
    - `->` for directed edges
    - Supports node, edge, and graph attributes
    - Supports parallel edges
    """
    # Local import to avoid circular dependencies
    from ...primitives.m_multigraph.io import to_phylozoo_dot as mm_to_phylozoo_dot
    
    # Use the underlying MixedMultiGraph's to_phylozoo_dot function
    return mm_to_phylozoo_dot(sd_network._graph, **kwargs)


def from_phylozoo_dot(pzdot_string: str, **kwargs: Any) -> 'SemiDirectedPhyNetwork':
    """
    Parse a PhyloZoo-DOT format string and create a SemiDirectedPhyNetwork.
    
    This function parses the PhyloZoo-DOT string into a MixedMultiGraph,
    then converts it to a SemiDirectedPhyNetwork.
    
    Parameters
    ----------
    pzdot_string : str
        PhyloZoo-DOT format string containing graph data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    SemiDirectedPhyNetwork
        Parsed semi-directed phylogenetic network.
    
    Raises
    ------
    ValueError
        If the PhyloZoo-DOT string is malformed or cannot be parsed, or if the resulting
        network structure is invalid for SemiDirectedPhyNetwork.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork.io import from_phylozoo_dot
    >>> 
    >>> pzdot_str = '''graph {
    ...     1 [label="A"];
    ...     2 [label="B"];
    ...     3;
    ...     1 -- 3;
    ...     2 -- 3;
    ... }'''
    >>> 
    >>> net = from_phylozoo_dot(pzdot_str)
    >>> net.number_of_nodes()
    3
    
    Notes
    -----
    The conversion process:
    1. Parse PhyloZoo-DOT string to MixedMultiGraph using from_phylozoo_dot
    2. Convert MixedMultiGraph to SemiDirectedPhyNetwork using sdnetwork_from_graph
    """
    # Local imports to avoid circular dependencies
    from .conversions import sdnetwork_from_graph
    from ...primitives.m_multigraph.io import from_phylozoo_dot as mm_from_phylozoo_dot
    
    # Parse PhyloZoo-DOT string to MixedMultiGraph
    mm_graph = mm_from_phylozoo_dot(pzdot_string, **kwargs)
    
    # Clean up labels: convert integer labels to strings, or remove labels that are just the node ID
    for node in mm_graph.nodes():
        # Check both directed and undirected graphs for node attributes
        if node in mm_graph._directed.nodes:
            node_data = mm_graph._directed.nodes[node]
            if 'label' in node_data:
                label = node_data['label']
                # If label is an integer and equals the node ID, remove it (auto-labeling will handle it)
                if isinstance(label, int) and label == node:
                    del node_data['label']
                # If label is an integer but different from node ID, convert to string
                elif isinstance(label, int):
                    node_data['label'] = str(label)
        
        if node in mm_graph._undirected.nodes:
            node_data = mm_graph._undirected.nodes[node]
            if 'label' in node_data:
                label = node_data['label']
                # If label is an integer and equals the node ID, remove it (auto-labeling will handle it)
                if isinstance(label, int) and label == node:
                    del node_data['label']
                # If label is an integer but different from node ID, convert to string
                elif isinstance(label, int):
                    node_data['label'] = str(label)
        
        if node in mm_graph._combined.nodes:
            node_data = mm_graph._combined.nodes[node]
            if 'label' in node_data:
                label = node_data['label']
                # If label is an integer and equals the node ID, remove it (auto-labeling will handle it)
                if isinstance(label, int) and label == node:
                    del node_data['label']
                # If label is an integer but different from node ID, convert to string
                elif isinstance(label, int):
                    node_data['label'] = str(label)
    
    # Convert MixedMultiGraph to SemiDirectedPhyNetwork
    return sdnetwork_from_graph(mm_graph, network_type='semi-directed')


# Register format handlers with FormatRegistry
# Use local import to avoid circular dependencies
def _register_formats() -> None:
    """Register format handlers for SemiDirectedPhyNetwork."""
    from .sd_phynetwork import SemiDirectedPhyNetwork
    
    FormatRegistry.register(
        SemiDirectedPhyNetwork, 'newick',
        reader=from_newick,
        writer=to_newick,
        extensions=['.nwk', '.newick', '.enewick', '.eNewick', '.enw'],
        default=True  # Newick is the default format
    )
    
    FormatRegistry.register(
        SemiDirectedPhyNetwork, 'phylozoo-dot',
        reader=from_phylozoo_dot,
        writer=to_phylozoo_dot,
        extensions=['.pzdot'],
        default=False
    )

# Register formats when module is imported
_register_formats()
