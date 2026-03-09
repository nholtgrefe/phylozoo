"""
Base module for semi-directed level-k generators.

This module provides the SemiDirectedGenerator class for representing level-k generators
of semi-directed phylogenetic networks. Generators are minimal biconnected components
that represent the core structure of level-k networks. These are not networks
themselves, but simplified structures used to build networks.

Unlike DirectedGenerator, SemiDirectedGenerator allows both directed and undirected edges.
"""

from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import (
    PhyloZooGeneratorStructureError, 
    PhyloZooValueError,
    PhyloZooGeneratorError,
)

from ....primitives.m_multigraph import MixedMultiGraph
from ....primitives.m_multigraph.features import (
    has_self_loops,
    bi_edge_connected_components,
    source_components,
)
from ....primitives.m_multigraph.transformations import orient_away_from_vertex
from ....primitives.d_multigraph import DirectedMultiGraph
from .....utils.validation import validation_aware, no_validation
from ...dnetwork.generator.side import Side, HybridSide, DirEdgeSide, IsolatedNodeSide
from ...dnetwork.generator.base import DirectedGenerator
from .side import BidirectedEdgeSide, UndirEdgeSide

if TYPE_CHECKING:
    from ...sdnetwork.base import MixedPhyNetwork

T = TypeVar('T')


@validation_aware(allowed=["validate", "_validate_*"], default=["validate"])
class SemiDirectedGenerator:
    """
    A level-k generator for semi-directed phylogenetic networks.
    
    A generator is a biconnected component that represents the
    core structure of a level-k network. Generators use MixedMultiGraph
    directly and have their own validation rules. These are not networks
    themselves, but simplified structures used to build networks.
    
    Unlike DirectedGenerator, SemiDirectedGenerator allows both directed
    and undirected edges.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The underlying graph structure of the generator.         Should be biconnected.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> # Create a generator with both directed and undirected edges
    >>> gen_graph = MixedMultiGraph(
    ...     directed_edges=[(0, 1), (0, 1)],  # Parallel directed edges
    ...     undirected_edges=[(1, 2)]  # Undirected edge
    ... )
    >>> generator = SemiDirectedGenerator(gen_graph)
    >>> generator.level
    1
    >>> generator.hybrid_nodes
    {1}
    
    Attributes
    ----------
    _graph : MixedMultiGraph[T]
        Internal graph structure using MixedMultiGraph.
        **Warning:** Do not modify directly.
    """
    
    def __init__(self, graph: MixedMultiGraph[T]) -> None:
        """
        Initialize a generator from a graph.
        
        Parameters
        ----------
        graph : MixedMultiGraph
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
    def graph(self) -> MixedMultiGraph[T]:
        """Get the underlying graph structure."""
        return self._graph
    
    @cached_property
    def hybrid_nodes(self) -> set[T]:
        """
        Get all hybrid nodes of this generator.

        A hybrid node is a node with in-degree >= 2 (from directed edges).
        Level-1 generators are one node with one undirected self-loop (bidirected edge);
        that node is the single hybrid.

        Returns
        -------
        set[T]
            Set of all hybrid node identifiers.

        Examples
        --------
        >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
        >>> gen_graph = MixedMultiGraph(directed_edges=[(0, 1), (0, 1)])
        >>> generator = SemiDirectedGenerator(gen_graph)
        >>> hybrid_nodes = generator.hybrid_nodes
        >>> 1 in hybrid_nodes
        True
        """
        # Level-1: one node with one undirected self-loop (bidirected edge)
        if self._graph.number_of_nodes() == 1 and self._graph.number_of_edges() == 1:
            node = list(self._graph.nodes())[0]
            if self._graph._undirected.number_of_edges(node, node) == 1:
                return {node}
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
        >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
        >>> gen_graph = MixedMultiGraph(directed_edges=[(0, 1), (0, 1)])
        >>> generator = SemiDirectedGenerator(gen_graph)
        >>> generator.level
        1
        """
        return len(self.hybrid_nodes)
    
    def validate(self) -> None:
        """
        Validate the generator structure.
        
        A semi-directed generator is valid if it can be rooted on an edge as a directed generator.
        This is done by:

        1. Finding source components
        2. Finding an edge in the source component
        3. Subdividing that edge
        4. Orienting away from the subdivision vertex
        5. Checking if the result is a valid DirectedGenerator
        
        Level-0: single node, no edges. Level-1: single node, one undirected
        self-loop (bidirected edge). Level >= 2: structural constraints and
        rootability as above.

        Raises
        ------
        PhyloZooGeneratorStructureError
            If any validation constraint is violated.
        """
        if self._graph.number_of_nodes() == 0:
            raise PhyloZooGeneratorStructureError(
                "Generator cannot be empty (must have at least one node)."
            )

        # Level-0: single node, no edges. Level-1: single node, one undirected self-loop.
        if self._graph.number_of_nodes() == 1:
            node = list(self._graph.nodes())[0]
            num_edges = self._graph.number_of_edges()
            if num_edges == 0:
                return  # Level-0
            if num_edges == 1:
                undir_self = list(
                    self._graph._undirected.edges(node, node, keys=True)
                )
                if len(undir_self) == 1:
                    return  # Level-1 (bidirected edge)
                raise PhyloZooGeneratorStructureError(
                    "Single node generator with one edge must have one "
                    "undirected self-loop (level-1)."
                )
            raise PhyloZooGeneratorStructureError(
                f"Single node generator must have 0 or 1 edge, got {num_edges}."
            )

        # Level >= 2: validate structural constraints and rootability
        self._validate_structural_constraints()
        self._validate_rootability()
    
    def _validate_structural_constraints(self) -> None:
        """
        Validate structural constraints (bi-edge connected, no self-loops, acyclic).
        
        Raises
        ------
        PhyloZooGeneratorStructureError
            If bi-edge connected, self-loop, or acyclicity constraints are violated.
        """
        # 1. Check that generator is bi-edge connected
        bi_edge_comps = list(bi_edge_connected_components(self._graph))
        if len(bi_edge_comps) > 1:
            raise PhyloZooGeneratorStructureError(
                f"Generator graph must be a single bi-edge connected component (blob), "
                f"but found {len(bi_edge_comps)} bi-edge connected components"
            )
        
        # 2. Disallow self-loops (except the special case handled in validate())
        if has_self_loops(self._graph):
            raise PhyloZooGeneratorStructureError("Self-loops are not allowed in SemiDirectedGenerator.")
        
        # 3. Check for directed cycles (must be acyclic in directed part)
        # Create a directed subgraph for cycle checking
        if self._graph._directed.number_of_edges() > 0:
            if not nx.is_directed_acyclic_graph(self._graph._directed):
                cycles = list(nx.simple_cycles(self._graph._directed))
                raise PhyloZooGeneratorStructureError(
                    f"Generator contains directed cycles. Found {len(cycles)} cycle(s). "
                    f"First cycle: {cycles[0] if cycles else 'unknown'}"
                )
    
    def _validate_rootability(self) -> None:
        """
        Validate that the generator can be rooted on an edge as a directed generator.
        
        This is done by:

        1. Finding source components
        2. Finding an edge in the source component
        3. Subdividing that edge
        4. Orienting away from the subdivision vertex
        5. Checking if the result is a valid DirectedGenerator
        
        Raises
        ------
        PhyloZooGeneratorStructureError
            If the generator cannot be rooted on an edge to form a valid DirectedGenerator.
        """
        # Step 1: Find source components
        components = source_components(self._graph)
        if len(components) != 1:
            raise PhyloZooGeneratorStructureError(
                f"Semi-directed generator must have exactly one source component; found {len(components)}"
            )
        nodes_in_component, undirected_edges_in_comp, outgoing_edges = components[0]
        if not nodes_in_component:
            raise PhyloZooGeneratorStructureError("Source component is empty")
        
        # Step 2: Find an edge in the source component.
        # We can only root on undirected edges in the source component or on
        # outgoing directed edges (from the component to outside).
        edge_to_subdivide = None
        edge_type: str | None = None

        # Prefer undirected edges inside the component
        if undirected_edges_in_comp:
            u, v, key = undirected_edges_in_comp[0]
            edge_to_subdivide = (u, v, key)
            edge_type = "undirected"
        # Otherwise, fall back to an outgoing directed edge
        elif outgoing_edges:
            u, v, key = outgoing_edges[0]
            edge_to_subdivide = (u, v, key)
            edge_type = "directed"
        
        if edge_to_subdivide is None:
            raise PhyloZooGeneratorStructureError(
                "No edge found in source component to subdivide for rooting"
            )
        
        u, v, key = edge_to_subdivide

        # Step 3: Subdivide the edge
        # Copy the whole graph, then remove the edge and insert the subdivision node.
        graph_copy = self._graph.copy()

        # Remove the edge to be subdivided.
        # For undirected edges, normalize (u, v) first so we match storage order.
        if edge_type == "undirected":
            norm_u, norm_v, _ = MixedMultiGraph.normalize_undirected_edge(u, v, key=key)
            graph_copy.remove_edge(norm_u, norm_v, key=key)
        else:
            graph_copy.remove_edge(u, v, key=key)
        
        # Add subdivision vertex using MixedMultiGraph's node-id generator
        subdiv_node = next(graph_copy.generate_node_ids(1))
        graph_copy.add_node(subdiv_node)
        
        # Add subdivision edges
        graph_copy.add_directed_edge(subdiv_node, u)
        graph_copy.add_directed_edge(subdiv_node, v)
        
        # Step 4: Orient away from subdivision vertex
        try:
            oriented_dm = orient_away_from_vertex(graph_copy, subdiv_node)
        except PhyloZooValueError as e:
            raise PhyloZooGeneratorStructureError(
                f"Failed to orient generator away from subdivision vertex: {e}"
            )
        
        # Step 5: Check if the result is a valid DirectedGenerator
        try:
            _ = DirectedGenerator(oriented_dm)
            # If we get here, it's valid
        except PhyloZooGeneratorError as e:
            raise PhyloZooGeneratorStructureError(
                f"Semi-directed generator cannot be rooted on edge to form valid "
                f"DirectedGenerator: {e}"
            )
    
    @cached_property
    def parallel_directed_edge_sides(self) -> list[tuple[DirEdgeSide, ...]]:
        """
        Get tuples of parallel directed edge sides.
        
        Returns groups of DirEdgeSide objects that represent parallel directed edges
        (multiple directed edges between the same pair of nodes).
        
        Returns
        -------
        list[tuple[DirEdgeSide, ...]]
            List of tuples, where each tuple contains DirEdgeSide objects
            representing parallel directed edges between the same pair of nodes.
        """
        # Group directed edges by (u, v) pair
        edge_groups: dict[tuple[T, T], list[DirEdgeSide]] = {}
        for u, v, key in self._graph.directed_edges_iter(keys=True):
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
    def non_parallel_directed_edge_sides(self) -> list[DirEdgeSide]:
        """
        Get all non-parallel directed edge sides of this generator.

        Returns all directed edges that are not part of parallel edge groups.

        Returns
        -------
        list[DirEdgeSide]
            List of non-parallel directed edge sides.
        """
        all_edges: set[tuple[T, T, int]] = set(
            self._graph.directed_edges_iter(keys=True)
        )
        parallel_edges: set[tuple[T, T, int]] = set()
        for parallel_group in self.parallel_directed_edge_sides:
            for edge_side in parallel_group:
                parallel_edges.add((edge_side.u, edge_side.v, edge_side.key))
        non_parallel_edges = all_edges - parallel_edges
        return [DirEdgeSide(u=u, v=v, key=key) for u, v, key in non_parallel_edges]

    @cached_property
    def directed_edge_sides(self) -> list[DirEdgeSide]:
        """
        Get all directed edge sides of this generator.

        Returns all directed edges in the generator as DirEdgeSide objects, including
        both parallel and non-parallel edges.

        Returns
        -------
        list[DirEdgeSide]
            List of all directed edges as DirEdgeSide objects.
        """
        edge_sides_list: list[DirEdgeSide] = []
        for parallel_group in self.parallel_directed_edge_sides:
            edge_sides_list.extend(parallel_group)
        edge_sides_list.extend(self.non_parallel_directed_edge_sides)
        return edge_sides_list

    @cached_property
    def undirected_edge_sides(self) -> list[UndirEdgeSide]:
        """
        Get all undirected edge sides (excluding the level-1 bidirected self-loop).

        Level-0 has no edges. Level-1 has only the undirected self-loop, which
        is represented as BidirectedEdgeSide in edge_sides, not here.

        Returns
        -------
        list[UndirEdgeSide]
            List of undirected edges with u != v.
        """
        if self.level == 0:
            return []
        elif self.level == 1:
            return []  # only edge is the self-loop (bidirected), excluded
        return [
            UndirEdgeSide(u=u, v=v, key=key)
            for u, v, key in self._graph.undirected_edges_iter(keys=True)
            if u != v
        ]

    @cached_property
    def edge_sides(self) -> list[DirEdgeSide | UndirEdgeSide | BidirectedEdgeSide]:
        """
        Get all edge sides (directed, undirected, and bidirected).

        Level-0 has no edges. Level-1 has one undirected self-loop (bidirected edge).
        Level >= 2 has directed and undirected edge sides.

        Returns
        -------
        list[DirEdgeSide | UndirEdgeSide | BidirectedEdgeSide]
            List of all edge sides.
        """
        if self.level == 0:
            return []
        elif self.level == 1:
            node = list(self._graph.nodes())[0]
            edge = next(
                iter(self._graph._undirected.edges(node, node, keys=True))
            )
            key = edge[2]
            return [BidirectedEdgeSide(node=node, key=key)]
        return list(self.directed_edge_sides) + list(self.undirected_edge_sides)
    
    @cached_property
    def hybrid_sides(self) -> list[HybridSide]:
        """
        Get all hybrid sides of this generator.
        
        Returns all hybrid nodes (in-degree >= 2 from directed edges) with out-degree 0
        (from directed edges) as HybridSide objects.
        
        Returns
        -------
        list[HybridSide]
            List of all hybrid nodes with out-degree 0 as HybridSide objects.
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

        Level-0: single node → [IsolatedNodeSide(node)].
        Level-1: one node, one bidirected edge → hybrid_sides + edge_sides.
        Level >= 2: hybrid_sides + edge_sides.

        Returns
        -------
        list[Side]
            List of all sides.
        """
        if self.level == 0:
            node = list(self._graph.nodes())[0]
            return [IsolatedNodeSide(node)]
        elif self.level == 1:
            return list(self.hybrid_sides) + list(self.edge_sides)
        else:
            return list(self.hybrid_sides) + list(self.edge_sides)
    
    def __repr__(self) -> str:
        """String representation of the generator."""
        num_nodes = self._graph.number_of_nodes()
        num_edges = self._graph.number_of_edges()
        return (
            f"SemiDirectedGenerator(level={self.level}, "
            f"nodes={num_nodes}, edges={num_edges})"
        )

