"""
Base module for directed level-k generators.

This module provides the DirectedGenerator class for representing level-k generators
of directed phylogenetic networks. Generators are minimal biconnected components
that represent the core structure of level-k networks. These are not networks
themselves, but simplified structures used to build networks.

Based on Gambette, Berry, and Paul (2009): "The Structure of Level-k Phylogenetic Networks".
"""

from __future__ import annotations

import warnings
from functools import cached_property
from typing import TYPE_CHECKING, Iterator, TypeVar

import networkx as nx

from ....primitives.d_multigraph import DirectedMultiGraph
from ....primitives.d_multigraph.features import (
    has_self_loops,
    bi_edge_connected_components,
)
from ....primitives.d_multigraph.transformations import subgraph as dm_subgraph
from ...dnetwork.classifications import is_binary, has_parallel_edges
from ...dnetwork.features import blobs
from ...dnetwork._utils import _suppress_deg2_nodes
from .....utils.validation import validation_aware
from .side import Side, HybridSide, DirEdgeSide

if TYPE_CHECKING:
    from ...dnetwork import DirectedPhyNetwork

T = TypeVar('T')


@validation_aware(allowed=["validate", "_validate_*"], default=["validate"])
class DirectedGenerator:
    """
    A level-k generator for directed phylogenetic networks.
    
    A generator is a biconnected component that represents the
    core structure of a level-k network. Generators use DirectedMultiGraph
    directly and have their own validation rules. These are not networks
    themselves, but simplified structures used to build networks.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The underlying graph structure of the generator. Should be biconnected.
    
    Attributes
    ----------
    _graph : DirectedMultiGraph[T]
        Internal graph structure using DirectedMultiGraph.
        **Warning:** Do not modify directly.
    _sides : list[Side] | None
        Cached list of sides (attachment points). Computed lazily.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> # Create a level-1 generator (root with parallel edges to hybrid node)
    >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])  # Parallel edges
    >>> generator = DirectedGenerator(gen_graph)
    >>> generator.level
    1
    >>> generator.hybrid_nodes
    {4}
    >>> # Create a level-0 generator (single node)
    >>> gen_graph0 = DirectedMultiGraph()
    >>> gen_graph0.add_node(1)
    >>> generator0 = DirectedGenerator(gen_graph0)
    >>> generator0.level
    0
    """
    
    def __init__(self, graph: DirectedMultiGraph[T]) -> None:
        """
        Initialize a generator from a graph.
        
        Parameters
        ----------
        graph : DirectedMultiGraph
            The graph structure of the generator. Should be biconnected.
        
        Raises
        ------
        ValueError
            If validation fails.
        """
        self._graph = graph
        
        # Validate the generator structure
        self.validate()
    
    @property
    def graph(self) -> DirectedMultiGraph[T]:
        """Get the underlying graph structure."""
        return self._graph
    
    @cached_property
    def hybrid_nodes(self) -> set[T]:
        """
        Get all hybrid nodes of this generator.
        
        A hybrid node is a node with in-degree >= 2.
        
        Returns
        -------
        set[T]
            Set of all hybrid node identifiers.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> hybrid_nodes = generator.hybrid_nodes
        >>> 4 in hybrid_nodes
        True
        """
        return {
            v for v in self._graph.nodes()
            if self._graph.indegree(v) >= 2
        }
    
    @cached_property
    def level(self) -> int:
        """
        Get the level of this generator.
        
        The level is the number of hybrid nodes in the generator.
        
        Returns
        -------
        int
            The level of the generator.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> generator.level
        1
        """
        return len(self.hybrid_nodes)
    
    def validate(self) -> None:
        """
        Validate the generator structure.
        
        Checks that:
        - The generator is not empty
        - Single node generators have no self-loops
        - Structural constraints (biconnected, no self-loops, acyclic)
        - Degree constraints (single root, no degree-1 nodes, etc.)
        
        Raises
        ------
        ValueError
            If any validation constraint is violated.
        """
        # 1. Check if empty generator
        if self._graph.number_of_nodes() == 0:
            raise ValueError("Generator cannot be empty (must have at least one node).")
        
        # 2. Check if single node generator
        if self._graph.number_of_nodes() == 1:
            if has_self_loops(self._graph):
                raise ValueError("Self-loops are not allowed in DirectedGenerator.")
            return
        
        # 3. Validate structural constraints
        self._validate_structural_constraints()
        
        # 4. Validate degree constraints
        self._validate_degree_constraints()
    
    def _validate_structural_constraints(self) -> None:
        """
        Validate structural constraints (bi-edge connected, no self-loops, acyclic).
        
        Raises
        ------
        ValueError
            If bi-edge connected, self-loop, or acyclicity constraints are violated.
        """
        # 1. Check that generator is bi-edge connected
        bi_edge_comps = list(bi_edge_connected_components(self._graph))
        if len(bi_edge_comps) > 1:
            raise ValueError(
                f"Generator graph must be a single bi-edge connected component (blob), "
                f"but found {len(bi_edge_comps)} bi-edge connected components"
            )
        
        # 2. Disallow self-loops
        if has_self_loops(self._graph):
            raise ValueError("Self-loops are not allowed in DirectedGenerator.")
        
        # 3. Check for directed cycles (must be acyclic)
        if not nx.is_directed_acyclic_graph(self._graph._graph):
            cycles = list(nx.simple_cycles(self._graph._graph))
            raise ValueError(
                f"Generator contains directed cycles. Found {len(cycles)} cycle(s). "
                f"First cycle: {cycles[0] if cycles else 'unknown'}"
            )
    
    def _validate_degree_constraints(self) -> None:
        """
        Validate degree constraints (single root, node degree patterns, etc.).
        
        Raises
        ------
        ValueError
            If degree constraints are violated.
        """
        # 1. Check that there is a single root node
        # This will raise an error if there are 0 or multiple roots
        _ = self.root_node
        
        # 2. Check degree constraints in a single iteration
        for node in self._graph.nodes():
            in_degree = self._graph.indegree(node)
            out_degree = self._graph.outdegree(node)
            total_degree = self._graph.degree(node)
            
            # 1) Check for degree-1 nodes
            if total_degree == 1:
                raise ValueError(
                    f"Generator has a degree-1 node: {node}"
                )
            
            # 2) If degree-2: either in-degree=2 or out-degree=2
            if total_degree == 2:
                if in_degree != 2 and out_degree != 2:
                    raise ValueError(
                        f"Generator has a degree-2 node {node} that does not have "
                        f"in-degree=2 or out-degree=2 (in-degree={in_degree}, out-degree={out_degree})"
                    )
            
            # 3) If degree-3: either in-degree=2 and out-degree=1, or in-degree=1 and out-degree=2
            elif total_degree == 3:
                if not ((in_degree == 2 and out_degree == 1) or (in_degree == 1 and out_degree == 2)):
                    raise ValueError(
                        f"Generator has a degree-3 node {node} that does not have "
                        f"(in-degree=2, out-degree=1) or (in-degree=1, out-degree=2) "
                        f"(in-degree={in_degree}, out-degree={out_degree})"
                    )
            
            # 4) If degree > 3: raise error
            elif total_degree > 3:
                raise ValueError(
                    f"Generator has a node {node} with degree > 3: {total_degree} "
                    f"(in-degree={in_degree}, out-degree={out_degree})"
                )
    
    @cached_property
    def root_node(self) -> T:
        """
        Get the root node of the generator.
        
        The root node is the unique node with in-degree 0.
        
        Returns
        -------
        T
            The root node identifier.
        
        Raises
        ------
        ValueError
            If there is no root node or multiple root nodes.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> generator.root_node
        8
        """
        roots = [v for v in self._graph.nodes() if self._graph.indegree(v) == 0]
        if len(roots) == 0:
            raise ValueError("Generator has no root node")
        if len(roots) > 1:
            raise ValueError(f"Generator has multiple root nodes: {roots}")
        return roots[0]
    
    @cached_property
    def parallel_edge_sides(self) -> list[tuple[DirEdgeSide, ...]]:
        """
        Get tuples of parallel edge sides.
        
        Returns groups of DirEdgeSide objects that represent parallel edges
        (multiple edges between the same pair of nodes).
        
        Returns
        -------
        list[tuple[DirEdgeSide, ...]]
            List of tuples, where each tuple contains DirEdgeSide objects
            representing parallel edges between the same pair of nodes.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> parallel = generator.parallel_edge_sides
        >>> len(parallel) > 0
        True
        """
        # Group edges by (u, v) pair
        edge_groups: dict[tuple[T, T], list[DirEdgeSide]] = {}
        for u, v, key in self._graph.edges(keys=True):
            edge_key = (u, v)
            if edge_key not in edge_groups:
                edge_groups[edge_key] = []
            edge_groups[edge_key].append(DirEdgeSide(u=u, v=v, key=key))
        
        # Return only groups with more than one edge (parallel edges)
        parallel_groups: list[tuple[DirEdgeSide, ...]] = []
        for edge_group in edge_groups.values():
            if len(edge_group) > 1:
                parallel_groups.append(tuple(edge_group))
        
        return parallel_groups
    
    @cached_property
    def non_parallel_edge_sides(self) -> list[DirEdgeSide]:
        """
        Get all non-parallel edge sides of this generator.
        
        Returns all edges that are not part of parallel edge groups.
        
        Returns
        -------
        list[DirEdgeSide]
            List of non-parallel edge sides.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> # Level-1 generator with only parallel edges (no non-parallel)
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> non_parallel = generator.non_parallel_edge_sides
        >>> len(non_parallel) == 0
        True
        """
        # Get all edges
        all_edges: set[tuple[T, T, int]] = set(self._graph.edges(keys=True))
        
        # Get edges that are part of parallel groups
        parallel_edges: set[tuple[T, T, int]] = set()
        for parallel_group in self.parallel_edge_sides:
            for edge_side in parallel_group:
                parallel_edges.add((edge_side.u, edge_side.v, edge_side.key))
        
        # Return non-parallel edges
        non_parallel_edges = all_edges - parallel_edges
        return [DirEdgeSide(u=u, v=v, key=key) for u, v, key in non_parallel_edges]
    
    @cached_property
    def edge_sides(self) -> list[DirEdgeSide]:
        """
        Get all edge sides of this generator.
        
        Returns all edges in the generator as DirEdgeSide objects, including
        both parallel and non-parallel edges.
        
        Returns
        -------
        list[DirEdgeSide]
            List of all edges as DirEdgeSide objects.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> edge_sides = generator.edge_sides
        >>> len(edge_sides) > 0
        True
        """
        # Sum of parallel (unpacked) and non-parallel edge sides
        edge_sides_list: list[DirEdgeSide] = []
        
        # Add all parallel edge sides (unpacked)
        for parallel_group in self.parallel_edge_sides:
            edge_sides_list.extend(parallel_group)
        
        # Add non-parallel edge sides
        edge_sides_list.extend(self.non_parallel_edge_sides)
        
        return edge_sides_list
    
    @cached_property
    def hybrid_sides(self) -> list[HybridSide]:
        """
        Get all hybrid sides of this generator.
        
        Returns all hybrid nodes (in-degree >= 2) with out-degree 0 as HybridSide objects.
        
        Returns
        -------
        list[HybridSide]
            List of all hybrid nodes with out-degree 0 as HybridSide objects.
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> hybrid_sides = generator.hybrid_sides
        >>> len(hybrid_sides) > 0
        True
        """
        hybrid_sides_list: list[HybridSide] = []
        for node in self.hybrid_nodes:
            if self._graph.outdegree(node) == 0:
                hybrid_sides_list.append(HybridSide(node=node))
        return hybrid_sides_list
    
    @cached_property
    def sides(self) -> list[Side]:
        """
        Get all sides (attachment points) of this generator.
        
        Returns both edge sides and hybrid sides.
        
        Returns
        -------
        list[Side]
            List of all sides (both DirEdgeSide and HybridSide objects).
        
        Examples
        --------
        >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
        >>> generator = DirectedGenerator(gen_graph)
        >>> sides = generator.sides
        >>> len(sides) > 0
        True
        """
        return list(self.edge_sides) + list(self.hybrid_sides)
    
    def __repr__(self) -> str:
        """String representation of the generator."""
        num_nodes = self._graph.number_of_nodes()
        num_edges = self._graph.number_of_edges()
        return (
            f"DirectedGenerator(level={self.level}, "
            f"nodes={num_nodes}, edges={num_edges})"
        )


def generators_from_network(network: 'DirectedPhyNetwork') -> Iterator[DirectedGenerator]:
    """
    Extract generators from a directed phylogenetic network.
    
    This function extracts all internal blobs from the network and converts
    them into generators. The network must be binary and have no parallel edges.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to extract generators from.
    
    Yields
    ------
    DirectedGenerator
        A generator for each internal blob in the network.
    
    Raises
    ------
    ValueError
        If the network is not binary.
    
    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> # Level-1 network with a single hybrid node
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (8, 5), (8, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid node 4
    ...         (4, 7),          # Hybrid to tree node
    ...         (5, 1), (6, 2), (7, 3), (7, 9)  # Tree nodes to leaves
    ...     ],
    ...     nodes=[
    ...         (1, {'label': 'A'}), (2, {'label': 'B'}),
    ...         (3, {'label': 'C'}), (9, {'label': 'D'})
    ...     ]
    ... )
    >>> generators = list(generators_from_network(net))
    >>> len(generators) > 0
    True
    """
    # Check that network is binary
    if not is_binary(network):
        raise ValueError("Network must be binary to extract generators")
    
    # Check that network has no parallel edges
    if has_parallel_edges(network):
        warnings.warn(
            "Network has parallel edges. The original paper does not treat networks with "
            "parallel edges. Proceed with care when using this, not everything may work "
            "as expected.",
            UserWarning,
            stacklevel=2
        )
    
    # Get all internal blobs (excluding leaves)
    internal_blobs = blobs(network, trivial=True, leaves=False)
    
    # Process each blob
    for blob_nodes in internal_blobs:
        # Get subgraph for this blob
        blob_graph = dm_subgraph(network._graph, blob_nodes)
        
        # Suppress all degree-2 nodes
        _suppress_deg2_nodes(blob_graph, exclude_nodes=None)
        
        # Create and yield the generator
        yield DirectedGenerator(blob_graph)
