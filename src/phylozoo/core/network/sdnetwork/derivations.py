"""
Network derivations module.

This module provides functions to derive other data structures from semi-directed
and mixed phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

# TODO: Implement derivation functions here.
# These could include:
# - Split extraction
# - Quartet extraction
# - Distance calculations
# - Blobtree construction
# - Subnetwork extraction
# - Displayed tree extraction
from typing import Any

from . import MixedPhyNetwork
from .features import blobs
from .transformations import suppress_2_blobs
from .sd_phynetwork import SemiDirectedPhyNetwork
from ...primitives.m_multigraph.transformations import (
    identify_vertices as mm_identify_vertices,
    suppress_degree2_node as mm_suppress_degree2_node,
)
from .io import sdnetwork_from_mmgraph


def tree_of_blobs(network: MixedPhyNetwork) -> MixedPhyNetwork:
    """
    Create the tree-of-blobs by suppressing all 2-blobs and collapsing internal blobs.

    This function:
    1. Suppresses all 2-blobs using suppress_2_blobs
    2. Finds all internal blobs (blobs with more than 1 node, excluding leaves)
    3. For each internal blob, identifies all vertices with a single vertex
    4. Returns a new network representing the tree-of-blobs

    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network to transform.

    Returns
    -------
    MixedPhyNetwork
        A new network where each blob has been collapsed to a single vertex,
        forming a tree structure.

    Examples
    --------
    >>> # Create a semi-directed network with a hybrid
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> sdnet = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         (5, 4),
    ...         (6, 4)
    ...     ],
    ...     undirected_edges=[
    ...         (5, 3),
    ...         (5, 6),
    ...         (6, 7),
    ...         (4, 8),
    ...         (8, 1),
    ...         (8, 2)
    ...     ],
    ...     nodes=[
    ...         (3, {'label': 'C'}),
    ...         (7, {'label': 'D'}),
    ...         (1, {'label': 'A'}),
    ...         (2, {'label': 'B'})
    ...     ]
    ... )
    >>> tree_net = tree_of_blobs(sdnet)
    >>> is_tree(tree_net)
    True
    """
    # First suppress all 2-blobs using the existing function
    blob_network = suppress_2_blobs(network)

    # Find all internal blobs (more than 1 node, not containing only leaves)
    all_blobs = blobs(blob_network, trivial=False, leaves=False)

    # Work on a copy of the internal graph and collapse blobs
    working_graph = blob_network._graph.copy()
    for blob in all_blobs:
        if len(blob) > 1:
            # Sort to ensure deterministic behavior
            blob_sorted = sorted(blob)
            # Identify all vertices in the blob with the first vertex
            mm_identify_vertices(working_graph, blob_sorted)

    # Convert back to appropriate network type
    # Preserve the input type (SemiDirectedPhyNetwork or MixedPhyNetwork)
    network_type = 'semi-directed' if isinstance(network, SemiDirectedPhyNetwork) else 'mixed'
    return sdnetwork_from_mmgraph(working_graph, network_type=network_type)
