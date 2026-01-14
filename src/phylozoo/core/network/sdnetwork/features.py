"""
Network features module.

This module provides functions to extract and identify features of semi-directed
and mixed phylogenetic networks (e.g., blobs, omnians, etc.).
"""

from functools import lru_cache
from typing import Any, TypeVar, Union

from ...primitives.m_multigraph.features import (
    bi_edge_connected_components,
    cut_edges as graph_cut_edges,
    cut_vertices as graph_cut_vertices,
    source_components,
)
from .base import MixedPhyNetwork
from .sd_phynetwork import SemiDirectedPhyNetwork
from phylozoo.utils.exceptions import (
    PhyloZooValueError,
    PhyloZooNetworkStructureError,
)

T = TypeVar('T')

# Root location can be a node or an edge
RootLocation = Union[T, tuple[T, T, int]]


@lru_cache(maxsize=128)
def blobs(
    network: MixedPhyNetwork,
    trivial: bool = True,
    leaves: bool = True,
) -> list[set[T]]:
    """
    Get blobs of the network.
    
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
    
    Returns
    -------
    list[set[T]]
        List of sets of nodes forming each blob.
    
    Raises
    ------
    PhyloZooValueError
        If `trivial=False` and `leaves=True` (this combination is not possible
        since leaves are single-node components).
    
    Notes
    -----
    Blobs are computed as bi-edge connected components (2-edge-connected components).
    A bi-edge connected component is a maximal subgraph that remains connected
    after removing any single edge (i.e., has no cut-edges/bridges).
    
    Results are cached using LRU cache with maxsize=128.
    
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
    >>> len([b for b in blobs(net, trivial=False, leaves=False)])
    1
    >>> # Filtering: exclude blobs containing only leaves
    >>> len([b for b in blobs(net, leaves=False)])
    1
    """
    # Check for invalid parameter combination
    if not trivial and leaves:
        raise PhyloZooValueError(
            "Cannot have trivial=False and leaves=True: leaves are single-node "
            "components, so excluding trivial components would exclude all leaves."
        )
    
    leaves_set = network.leaves
    result: list[set[T]] = []
    
    # Process bi-edge connected components directly
    for blob in bi_edge_connected_components(network._graph):
        blob_set = set(blob)
        
        # Filter single-node components based on parameters
        if len(blob_set) == 1:
            node = next(iter(blob_set))
            # Skip if trivial=False and it's not a leaf
            if not trivial and node not in leaves_set:
                continue
            # Skip if leaves=False and it is a leaf
            if not leaves and node in leaves_set:
                continue
            result.append(blob_set)
        else:
            # Multi-node components: always include (cannot consist entirely of leaves)
            result.append(blob_set)
    
    return result


@lru_cache(maxsize=128)
def k_blobs(
    network: MixedPhyNetwork,
    k: int,
    trivial: bool = True,
    leaves: bool = True,
) -> list[set[T]]:
    """
    Get k-blobs of the network.
    
    A k-blob is a blob with exactly k edges incident to it. An incident edge
    is any edge (directed or undirected) that connects a node inside the blob
    to a node outside the blob. Parallel edges are counted separately.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network.
    k : int
        The number of edges that should be incident to each returned blob.
    trivial : bool, optional
        Whether to include trivial (single-node) blobs. By default True.
    leaves : bool, optional
        Whether to include blobs that contain only leaves. By default True.
    
    Returns
    -------
    list[set[T]]
        List of sets of nodes forming each k-blob.
    
    Raises
    ------
    PhyloZooValueError
        If `trivial=False` and `leaves=True` (this combination is not possible
        since leaves are single-node components).
    
    Notes
    -----
    This function identifies blobs and then filters
    them based on the number of incident edges. Both directed and undirected
    edges are counted. Parallel edges are counted separately, so if there are
    two parallel edges crossing the blob boundary, they count as two incident edges.
    
    Results are cached using LRU cache with maxsize=128.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> # Tree network: leaves are 1-blobs, internal nodes are 2-blobs or more
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> sorted([sorted(b) for b in k_blobs(net, k=1)])
    [[1], [2]]
    >>> sorted([sorted(b) for b in k_blobs(net, k=2)])
    [[3]]
    """
    # Check for invalid parameter combination
    if not trivial and leaves:
        raise PhyloZooValueError(
            "Cannot have trivial=False and leaves=True: leaves are single-node "
            "components, so excluding trivial components would exclude all leaves."
        )
    
    # Get the underlying graph for edge access
    graph = network._graph
    result: list[set[T]] = []
    
    # Iterate through blobs
    for blob in blobs(network, trivial=trivial, leaves=leaves):
        blob_set = set(blob)
        
        # Count edges incident to this blob
        # An incident edge has exactly one endpoint in the blob
        incident_edge_count = 0
        seen_undirected = set()  # Track undirected edges to avoid double-counting
        
        # Check all edges incident to nodes in the blob
        for node in blob_set:
            # Check directed incoming edges (parent edges)
            for u, v, key in graph.incident_parent_edges(node, keys=True):
                if u not in blob_set:  # Edge crosses blob boundary
                    incident_edge_count += 1
            
            # Check directed outgoing edges (child edges)
            for u, v, key in graph.incident_child_edges(node, keys=True):
                if v not in blob_set:  # Edge crosses blob boundary
                    incident_edge_count += 1
            
            # Check undirected edges
            for u, v, key in graph.incident_undirected_edges(node, keys=True):
                neighbor = v if u == node else u
                if neighbor not in blob_set:  # Edge crosses blob boundary
                    # Normalize edge to avoid double-counting undirected edges
                    edge = graph.normalize_undirected_edge(u, v, key)
                    if edge not in seen_undirected:
                        seen_undirected.add(edge)
                        incident_edge_count += 1
        
        # Add blob if it has exactly k incident edges
        if incident_edge_count == k:
            result.append(blob_set)
    
    return result


@lru_cache(maxsize=128)
def cut_edges(network: MixedPhyNetwork) -> set[tuple[T, T, int]]:
    """
    Find all cut-edges (bridges) in the network.
    
    A cut-edge is an edge whose removal increases the number of
    connected components. Results are cached per network instance.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network.
    
    Returns
    -------
    set[tuple[T, T, int]]
        Set of cut-edges as 3-tuples (u, v, key).
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
    ...     undirected_edges=[(4, 2), (5, 8), (6, 9), (5, 10), (6, 11)],
    ...     nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
    ... )
    >>> edges = cut_edges(net)
    >>> len(edges) > 0
    True
    
    Notes
    -----
    Results are cached using LRU cache with maxsize=128.
    """
    return graph_cut_edges(network._graph, keys=True, data=False)


@lru_cache(maxsize=128)
def cut_vertices(network: MixedPhyNetwork) -> set[T]:
    """
    Find all cut-vertices (articulation points) in the network.
    
    A cut-vertex is a vertex whose removal increases the number of
    connected components. Results are cached per network instance.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network.
    
    Returns
    -------
    set[T]
        Set of cut-vertices.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
    ...     undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
    ...     nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
    ... )
    >>> vertices = cut_vertices(net)
    >>> len(vertices) > 0
    True
    
    Notes
    -----
    Results are cached using LRU cache with maxsize=128.
    """
    return graph_cut_vertices(network._graph, data=False)


def omnians(network: MixedPhyNetwork) -> set[T]:
    """Stub for omnians function."""
    return set()


def root_locations(
    network: SemiDirectedPhyNetwork,
) -> tuple[list[T], list[tuple[T, T, int]], list[tuple[T, T, int]]]:
    """
    Get all valid root locations for a semi-directed network.
    
    Returns three separate lists:
    1. Node root locations: non-leaf nodes in the source component
    2. Undirected edge root locations: undirected edges in the source component
    3. Directed edge root locations: directed edges exiting the source component
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    
    Raises
    ------
    PhyloZooNetworkStructureError
        If the network has multiple source components or an empty source component.
    
    Returns
    -------
    tuple[list[T], list[tuple[T, T, int]], list[tuple[T, T, int]]]
        A tuple containing:
        - List of node root locations (non-leaf nodes in source component)
        - List of undirected edge root locations as (u, v, key) tuples
        - List of directed edge root locations as (u, v, key) tuples (outgoing from source component)
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> node_locs, undir_locs, dir_locs = root_locations(net)
    >>> 3 in node_locs  # Node 3 is a valid root location
    True
    >>> (3, 1, 0) in undir_locs  # Edge (3, 1) is also a valid root location
    True
    """
    # Find source components
    components = source_components(network._graph)
    if len(components) != 1:
        raise PhyloZooNetworkStructureError(
            f"Semi-directed network must have exactly one source component; found {len(components)}"
        )
    
    nodes_in_component, undirected_edges_in_comp, outgoing_edges = components[0]
    if not nodes_in_component:
        raise PhyloZooNetworkStructureError("Source component is empty")
    
    internal_nodes = network.internal_nodes
    
    # 1. Node root locations: all non-leaf nodes in the source component
    node_locations: list[T] = [
        node for node in nodes_in_component if node in internal_nodes
    ]
    
    # 2. Undirected edge root locations: all undirected edges in the source component
    undirected_edge_locations: list[tuple[T, T, int]] = list(undirected_edges_in_comp)
    
    # 3. Directed edge root locations: exactly the outgoing edges from the source component
    directed_edge_locations: list[tuple[T, T, int]] = list(outgoing_edges)
    
    return (node_locations, undirected_edge_locations, directed_edge_locations)

