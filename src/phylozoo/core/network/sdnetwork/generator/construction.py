"""
Construction module for semi-directed level-k generators.

This module provides functions for constructing semi-directed generators
from directed generators.
"""

from __future__ import annotations

from typing import Any

from ....primitives.m_multigraph import MixedMultiGraph
from ....primitives.m_multigraph.isomorphism import is_isomorphic, _get_graph_invariant
from ....primitives.m_multigraph.transformations import suppress_degree2_node
from ...dnetwork.generator.base import DirectedGenerator
from ...dnetwork.generator.construction import all_level_k_generators as all_level_k_dgenerators
from .base import SemiDirectedGenerator
from phylozoo.utils.exceptions import PhyloZooValueError


def dgenerator_to_sdgenerator(d_generator: DirectedGenerator) -> SemiDirectedGenerator:
    """
    Convert a DirectedGenerator to a SemiDirectedGenerator.
    
    This function semi-directs a directed generator by:

    1. Converting all edges to undirected, except those entering a hybrid node
       (a node with in-degree >= 2).
    2. Suppressing the degree-2 root node.
    
    Parameters
    ----------
    d_generator : DirectedGenerator
        The directed generator to convert.
    
    Returns
    -------
    SemiDirectedGenerator
        A new semi-directed generator with non-hybrid edges undirected
        and the degree-2 root suppressed.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.core.network.dnetwork.generator.base import DirectedGenerator
    >>> # Create a level-1 directed generator (root with parallel edges to hybrid)
    >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
    >>> d_gen = DirectedGenerator(gen_graph)
    >>> sd_gen = dgenerator_to_sdgenerator(d_gen)
    >>> sd_gen.level
    1
    >>> # Level-1 is represented as one node with one undirected self-loop (bidirected edge)
    >>> sd_gen.graph.number_of_nodes()
    1
    >>> sd_gen.graph.number_of_edges()
    1
    >>> list(sd_gen.graph.undirected_edges_iter(keys=True))
    [(4, 4, 0)]
    """
    # Level-0: single node, no edges
    if d_generator.level == 0:
        mixed_graph = MixedMultiGraph()
        mixed_graph.add_node(d_generator.root_node)
        return SemiDirectedGenerator(mixed_graph)

    # Level-1 directed generator: two nodes, two parallel edges (root -> hybrid).
    # Represent as one node with one undirected self-loop (bidirected edge).
    if d_generator.level == 1:
        (hybrid_node,) = d_generator.hybrid_nodes
        return SemiDirectedGenerator(
            MixedMultiGraph(undirected_edges=[(hybrid_node, hybrid_node)])
        )

    # Get hybrid nodes (nodes with in-degree >= 2)
    hybrid_nodes = d_generator.hybrid_nodes

    # Collect all edges and their data, separating into directed and undirected
    directed_edges: list[dict[str, Any]] = []
    undirected_edges: list[dict[str, Any]] = []
    
    for u, v, key, data in d_generator.graph.edges(keys=True, data=True):
        edge_dict: dict[str, Any] = {"u": u, "v": v}
        if key != 0:
            edge_dict["key"] = key
        if data:
            edge_dict.update(data)
        
        # If v is a hybrid node, keep edge as directed; otherwise, make it undirected
        if v in hybrid_nodes:
            directed_edges.append(edge_dict)
        else:
            undirected_edges.append(edge_dict)
    
    # Create MixedMultiGraph with the separated edges
    mixed_graph = MixedMultiGraph(
        directed_edges=directed_edges,
        undirected_edges=undirected_edges
    )
    
    # Suppress degree-2 root node
    # The root node in a generator has in-degree 0 and typically out-degree 2
    # After conversion, if it's degree-2, it should be suppressed
    root_node = d_generator.root_node
    if root_node in mixed_graph.nodes() and mixed_graph.degree(root_node) == 2:
        suppress_degree2_node(mixed_graph, root_node, merged_attrs=None)

    # Create and return SemiDirectedGenerator
    return SemiDirectedGenerator(mixed_graph)


def all_level_k_generators(k: int) -> set[SemiDirectedGenerator]:
    """
    Generate all (strict) level-k semi-directed generators.
    
    This function constructs all (strict) level-k semi-directed generators by:

    1. Getting all level-k directed generators
    2. Semi-directing each one
    3. Filtering out isomorphic generators
    
    Parameters
    ----------
    k : int
        The level of generators to generate.
    
    Returns
    -------
    set[SemiDirectedGenerator]
        Set of all level-k semi-directed generators (up to isomorphism).
    
    Raises
    ------
    PhyloZooValueError
        If level is negative.
        
    Notes
    -----
    Validation is deferred during construction for performance optimization.
    The algorithm provably produces valid generators, so validation is skipped
    by default. To validate a generator, call ``gen.validate()`` explicitly.
    
    Examples
    --------
    >>> # Get all level-1 semi-directed generators
    >>> level1 = all_level_k_generators(1)
    >>> len(level1) == 1
    True
    >>> # Get all level-2 semi-directed generators
    >>> level2 = all_level_k_generators(2)
    >>> len(level2) == 2
    True
    """
    if k < 0:
        raise PhyloZooValueError("Level must be non-negative")
    
    # Get all level-k directed generators
    d_generators = all_level_k_dgenerators(k)
    
    # Semi-direct all of them
    result: list[SemiDirectedGenerator] = []
    # Dictionary mapping invariant -> list of generators with that invariant
    invariant_groups: dict[tuple[int, int, tuple[int, ...], tuple[int, ...], tuple[int, ...], tuple[int, ...]], list[SemiDirectedGenerator]] = {}
    
    for d_gen in d_generators:
        sd_gen = dgenerator_to_sdgenerator(d_gen)
        
        # Compute invariants for fast filtering
        invariant = _get_graph_invariant(sd_gen.graph)
        
        # Only check isomorphism against generators with matching invariants
        candidate_group = invariant_groups.get(invariant)
        is_duplicate = False
        
        if candidate_group is not None:
            # Only check against candidates with matching invariants
            for existing_gen in candidate_group:
                if is_isomorphic(sd_gen.graph, existing_gen.graph):
                    is_duplicate = True
                    break
        
        # Only add if not isomorphic to any existing generator
        if not is_duplicate:
            result.append(sd_gen)
            # Add to invariant group for future comparisons
            if candidate_group is None:
                invariant_groups[invariant] = [sd_gen]
            else:
                candidate_group.append(sd_gen)
    
    return set(result)

