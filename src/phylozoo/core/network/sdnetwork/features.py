"""
Network features module.

This module provides functions to extract and identify features of semi-directed
and mixed phylogenetic networks (e.g., blobs, omnians, etc.).
"""

from typing import Iterator, TypeVar

from ...primitives.m_multigraph.features import biconnected_components
from .base import MixedPhyNetwork

T = TypeVar('T')


def blobs(
    network: MixedPhyNetwork,
    trivial: bool = True,
    leaves: bool = True,
) -> Iterator[set[T]]:
    """
    Yield blobs of the network.
    
    A blob is a maximal subgraph without any cut-edges. This function provides
    filtering options to control which blobs are returned.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network.
    trivial : bool, optional
        Whether to include trivial (single-node) blobs. By default True.
    leaves : bool, optional
        Whether to include blobs that contain only leaves. By default True.
    
    Yields
    ------
    set[T]
        Sets of nodes forming each blob.
    
    Raises
    ------
    ValueError
        If `trivial=False` and `leaves=True` (this combination is not possible
        since leaves are single-node components).
    
    Notes
    -----
    Blobs are computed as:
    - All biconnected components with at least 3 nodes
    - Biconnected components with exactly 2 nodes that have parallel edges
    - All remaining nodes (trivial blobs, including leaves)
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> # Network with hybrid node creating a non-trivial blob
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         {'u': 6, 'v': 5, 'gamma': 0.6},  # Hybrid edge to node 5
    ...         {'u': 7, 'v': 5, 'gamma': 0.4},  # Hybrid edge to node 5
    ...     ],
    ...     undirected_edges=[
    ...         (5, 1),  # Hybrid node 5: in-degree 2, total degree 3 (2+1)
    ...         (6, 2), (6, 3), (6, 7),  # Tree node 6 connects to leaves and node 7
    ...         (7, 8), (7, 9),  # Tree node 7 connects to leaves
    ...     ],
    ...     nodes=[
    ...         (1, {'label': 'A'}),
    ...         (2, {'label': 'B'}),
    ...         (3, {'label': 'C'}),
    ...         (8, {'label': 'D'}),
    ...         (9, {'label': 'E'}),
    ...     ]
    ... )
    >>> sorted([sorted(b) for b in blobs(net)])
    [[1], [2], [3], [5, 6, 7], [8], [9]]
    >>> # Filtering: exclude trivial (single-node) blobs
    >>> len(list(blobs(net, trivial=False, leaves=False)))
    1
    >>> # Filtering: exclude blobs containing only leaves
    >>> len(list(blobs(net, leaves=False)))
    1
    """
    # Check for invalid parameter combination
    if not trivial and leaves:
        raise ValueError(
            "Cannot have trivial=False and leaves=True: leaves are single-node "
            "components, so excluding trivial components would exclude all leaves."
        )
    
    leaves_set = network.leaves
    graph = network._graph
    
    # Collect nodes that are in non-trivial blobs
    visited: set[T] = set()
    
    # Process biconnected components
    for blob in biconnected_components(graph):
        blob_set = set(blob)
        
        # Keep blobs with >= 3 nodes
        if len(blob_set) >= 3:
            visited.update(blob_set)
            # Filter blobs containing only leaves
            if not leaves and blob_set.issubset(leaves_set):
                continue
            yield blob_set
            continue
        
        # Keep 2-node blobs only if they have parallel edges
        if len(blob_set) == 2:
            u, v = list(blob_set)
            # Check for parallel edges in both directed and undirected graphs
            num_edges = 0
            # Count directed edges (both directions)
            if u in graph._directed and v in graph._directed[u]:
                num_edges += len(graph._directed[u][v])
            if v in graph._directed and u in graph._directed[v]:
                num_edges += len(graph._directed[v][u])
            # Count undirected edges
            if u in graph._undirected and v in graph._undirected[u]:
                num_edges += len(graph._undirected[u][v])
            
            if num_edges > 1:  # Parallel edges exist
                visited.update(blob_set)
                yield blob_set
    
    # Add remaining nodes as trivial blobs
    for node in graph.nodes():
        if node not in visited:
            if node not in leaves_set and trivial:
                yield {node}
            elif node in leaves_set and leaves:
                yield {node}

