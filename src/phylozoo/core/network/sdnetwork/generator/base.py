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

from ....primitives.m_multigraph import MixedMultiGraph
from ....primitives.m_multigraph.features import (
    has_self_loops,
    bi_edge_connected_components,
    source_components,
)
from ....primitives.m_multigraph.transformations import orient_away_from_vertex
from ....primitives.d_multigraph import DirectedMultiGraph
from .....utils.validation import validation_aware, no_validation
from ...dnetwork.generator.side import Side, HybridSide, DirEdgeSide
from ...dnetwork.generator.base import DirectedGenerator
from .side import UndirEdgeSide

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
        The underlying graph structure of the generator. Should be biconnected.
    
    Attributes
    ----------
    _graph : MixedMultiGraph[T]
        Internal graph structure using MixedMultiGraph.
        **Warning:** Do not modify directly.
    
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
        
        Special case: A level-1 generator can be a single node with two directed self-loops
        (representing one bidirected self-loop), which should be valid and have one hybrid.
        
        Raises
        ------
        ValueError
            If any validation constraint is violated.
        """
        # 1. Check if empty generator
        if self._graph.number_of_nodes() == 0:
            raise ValueError("Generator cannot be empty (must have at least one node).")
        
        # 2. Special case: single node with two directed self-loops (level-1 generator)
        if self._graph.number_of_nodes() == 1:
            node = list(self._graph.nodes())[0]
            # Check for two directed self-loops
            self_loop_count = self._graph._directed.number_of_edges(node, node)
            if self_loop_count == 2:
                # This is a valid level-1 generator
                # Check that it has exactly one hybrid (the node itself)
                if len(self.hybrid_nodes) == 1 and node in self.hybrid_nodes:
                    return  # Valid - skip other validations
                else:
                    raise ValueError(
                        "Single node with two directed self-loops must have exactly one hybrid node"
                    )
            elif self_loop_count > 0:
                raise ValueError(
                    f"Single node generator has {self_loop_count} self-loops, "
                    "but only 2 directed self-loops are allowed for level-1 generators"
                )
            # No self-loops - this is a level-0 generator (single node, no edges)
            # But we still need to check if it has any edges
            if self._graph.number_of_edges() > 0:
                raise ValueError(
                    "Single node generator with edges must have exactly 2 directed self-loops"
                )
            return  # Valid level-0 generator
        
        # 3. Validate structural constraints
        self._validate_structural_constraints()
        
        # 4. Validate that it can be rooted on an edge
        self._validate_can_be_rooted_on_edge()
    
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
        
        # 2. Disallow self-loops (except the special case handled in validate())
        if has_self_loops(self._graph):
            raise ValueError("Self-loops are not allowed in SemiDirectedGenerator.")
        
        # 3. Check for directed cycles (must be acyclic in directed part)
        # Create a directed subgraph for cycle checking
        if self._graph._directed.number_of_edges() > 0:
            if not nx.is_directed_acyclic_graph(self._graph._directed):
                cycles = list(nx.simple_cycles(self._graph._directed))
                raise ValueError(
                    f"Generator contains directed cycles. Found {len(cycles)} cycle(s). "
                    f"First cycle: {cycles[0] if cycles else 'unknown'}"
                )
    
    def _validate_can_be_rooted_on_edge(self) -> None:
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
        ValueError
            If the generator cannot be rooted on an edge to form a valid DirectedGenerator.
        """
        # Step 1: Find source components
        components = source_components(self._graph)
        if len(components) != 1:
            raise ValueError(
                f"Semi-directed generator must have exactly one source component; found {len(components)}"
            )
        nodes_in_component, undirected_edges_in_comp, outgoing_edges = components[0]
        if not nodes_in_component:
            raise ValueError("Source component is empty")
        
        # Step 2: Find an edge in the source component
        # source_components returns (nodes, undirected_edges, outgoing_edges)
        # where undirected_edges are (u, v, key) tuples within the component
        # and outgoing_edges are (u, v, key) with u in component, v not in component
        edge_to_subdivide = None
        edge_type = None
        nodes_set = set(nodes_in_component)
        
        # Try undirected edges first (these are within the component)
        if undirected_edges_in_comp:
            # Take the first undirected edge
            u, v, key = undirected_edges_in_comp[0]
            edge_to_subdivide = (u, v, key)
            edge_type = 'undirected'
        else:
            # If no undirected edges, try directed edges within the component
            # Check all directed edges to find one where both endpoints are in the component
            for u, v, key in self._graph.directed_edges_iter(keys=True):
                if u in nodes_set and v in nodes_set:
                    edge_to_subdivide = (u, v, key)
                    edge_type = 'directed'
                    break
            
            # If still no edge found, try using an outgoing edge
            # (this can happen when the component is a single node with outgoing edges)
            if edge_to_subdivide is None and outgoing_edges:
                u, v, key = outgoing_edges[0]
                edge_to_subdivide = (u, v, key)
                edge_type = 'directed'
        
        if edge_to_subdivide is None:
            raise ValueError(
                "No edge found in source component to subdivide for rooting"
            )
        
        u, v, key = edge_to_subdivide
        
        # Step 3: Subdivide the edge
        # Create a copy of the graph
        graph_copy = MixedMultiGraph()
        
        # Add all nodes
        for node in self._graph.nodes():
            graph_copy.add_node(node)
        
        # Add all edges except the one to subdivide
        for edge_u, edge_v, edge_key in self._graph.directed_edges_iter(keys=True):
            if (edge_u, edge_v, edge_key) != (u, v, key):
                data = self._graph._directed[edge_u][edge_v][edge_key]
                graph_copy.add_directed_edge(edge_u, edge_v, key=edge_key, **data)
        
        for edge_u, edge_v, edge_key in self._graph.undirected_edges_iter(keys=True):
            if (edge_u, edge_v, edge_key) != (u, v, key) and (edge_v, edge_u, edge_key) != (u, v, key):
                data = self._graph._undirected[edge_u][edge_v][edge_key]
                graph_copy.add_undirected_edge(edge_u, edge_v, key=edge_key, **data)
        
        # Add subdivision vertex
        max_node = max(self._graph.nodes(), default=-1)
        subdiv_node = max_node + 1
        graph_copy.add_node(subdiv_node)
        
        # Add edges from subdivision: u -> subdiv_node -> v
        if edge_type == 'undirected':
            # For undirected edge, add u -> subdiv_node and subdiv_node -> v
            graph_copy.add_directed_edge(u, subdiv_node)
            graph_copy.add_directed_edge(subdiv_node, v)
        else:
            # For directed edge, preserve direction: u -> subdiv_node -> v
            graph_copy.add_directed_edge(u, subdiv_node)
            graph_copy.add_directed_edge(subdiv_node, v)
        
        # Step 4: Orient away from subdivision vertex
        try:
            oriented_dm = orient_away_from_vertex(graph_copy, subdiv_node)
        except ValueError as e:
            raise ValueError(
                f"Failed to orient generator away from subdivision vertex: {e}"
            )
        
        # Step 5: Check if the result is a valid DirectedGenerator
        with no_validation(classes=["DirectedGenerator"]):
            try:
                d_gen = DirectedGenerator(oriented_dm)
                # If we get here, it's valid
            except ValueError as e:
                raise ValueError(
                    f"Semi-directed generator cannot be rooted on edge to form valid "
                    f"DirectedGenerator: {e}"
                )
    
    def _validate_degree_constraints(self) -> None:
        """
        Validate degree constraints (node degree patterns, etc.).
        
        Raises
        ------
        ValueError
            If degree constraints are violated.
        
        Notes
        -----
        Validation is currently a stub and will be implemented in the future.
        """
        # TODO: Implement validation
        # Check degree constraints in a single iteration
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
        # Get all directed edges
        all_edges: set[tuple[T, T, int]] = set(self._graph.directed_edges_iter(keys=True))
        
        # Get edges that are part of parallel groups
        parallel_edges: set[tuple[T, T, int]] = set()
        for parallel_group in self.parallel_directed_edge_sides:
            for edge_side in parallel_group:
                parallel_edges.add((edge_side.u, edge_side.v, edge_side.key))
        
        # Return non-parallel edges
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
        
        # Add all parallel directed edge sides (unpacked)
        for parallel_group in self.parallel_directed_edge_sides:
            edge_sides_list.extend(parallel_group)
        
        # Add non-parallel directed edge sides
        edge_sides_list.extend(self.non_parallel_directed_edge_sides)
        
        return edge_sides_list
    
    @cached_property
    def undirected_edge_sides(self) -> list[UndirEdgeSide]:
        """
        Get all undirected edge sides of this generator.
        
        Semi-directed generators do not have parallel undirected edges, so this
        returns all undirected edges.
        
        Returns
        -------
        list[UndirEdgeSide]
            List of all undirected edges as UndirEdgeSide objects.
        """
        # Semi-directed generators have no parallel undirected edges
        # So all undirected edges are returned directly
        return [
            UndirEdgeSide(u=u, v=v, key=key)
            for u, v, key in self._graph.undirected_edges_iter(keys=True)
        ]
    
    @cached_property
    def edge_sides(self) -> list[DirEdgeSide | UndirEdgeSide]:
        """
        Get all edge sides (both directed and undirected) of this generator.
        
        Returns
        -------
        list[DirEdgeSide | UndirEdgeSide]
            List of all edge sides (both directed and undirected).
        """
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
        
        Returns both edge sides (directed and undirected) and hybrid sides.
        
        Returns
        -------
        list[Side]
            List of all sides (DirEdgeSide, UndirEdgeSide, and HybridSide objects).
        """
        return list(self.edge_sides) + list(self.hybrid_sides)
    
    def __repr__(self) -> str:
        """String representation of the generator."""
        num_nodes = self._graph.number_of_nodes()
        num_edges = self._graph.number_of_edges()
        return (
            f"SemiDirectedGenerator(level={self.level}, "
            f"nodes={num_nodes}, edges={num_edges})"
        )

