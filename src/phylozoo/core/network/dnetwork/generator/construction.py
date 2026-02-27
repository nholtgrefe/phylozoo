"""
Construction module for building level-k generators from level-(k-1) generators.

This module implements the R1 and R2 transformation rules from Gambette, Berry, and Paul (2009)
to construct all level-k generators from level-(k-1) generators.

Based on Gambette, Berry, and Paul (2009): "The Structure of Level-k Phylogenetic Networks".
"""

from __future__ import annotations

from typing import Any, Iterator

import networkx as nx

from ....primitives.d_multigraph import DirectedMultiGraph
from ....primitives.d_multigraph.isomorphism import is_isomorphic, _get_graph_invariant
from .....utils.validation import no_validation
from .base import DirectedGenerator
from .side import Side, DirEdgeSide, EdgeSide, HybridSide, NodeSide
from phylozoo.utils.exceptions import PhyloZooValueError


def _get_node_reachability_matrix(
    generator: DirectedGenerator,
) -> dict[tuple[Any, Any], bool]:
    """
    Compute reachability matrix for all pairs of nodes efficiently.
    
    Uses topological sorting and dynamic programming for DAGs.
    
    Parameters
    ----------
    generator : DirectedGenerator
        The generator to compute node reachability for.
    
    Returns
    -------
    dict[tuple[Any, Any], bool]
        Dictionary mapping (source_node, target_node) -> bool indicating
        if there's a directed path from source to target.
        Only stores True values; missing keys indicate False.
    
    Notes
    -----
    - Uses topological sort + dynamic programming for efficient DAG processing.
    - Each node is reachable from itself.
    - Only stores True values; use .get(key, False) to check reachability.
    """
    graph = generator.graph._graph  # NetworkX MultiDiGraph
    nodes = list(graph.nodes())
    
    # Compute node-to-node reachability using topological sort (DAG optimization)
    topo_order = list(nx.topological_sort(graph))
    
    # Use a set of reachable pairs for faster lookups, then convert to dict
    reachable_pairs: set[tuple[Any, Any]] = set()
    
    # Initialize: each node is reachable from itself
    for node in nodes:
        reachable_pairs.add((node, node))
    
    # Process in reverse topological order
    # For each node v, compute which nodes are reachable from v
    # by taking union of reachable nodes from all successors of v
    for v in reversed(topo_order):
        reachable_from_v: set[Any] = {v}
        
        # Add all nodes reachable from successors of v
        for successor in graph.successors(v):
            reachable_from_v.add(successor)
            # Add all nodes reachable from successor
            for node in nodes:
                if (successor, node) in reachable_pairs:
                    reachable_from_v.add(node)
        
        # Store reachability from v
        for target in reachable_from_v:
            reachable_pairs.add((v, target))
    
    # Convert to dict with True values only (missing keys indicate False)
    return {pair: True for pair in reachable_pairs}


def _apply_R1(
    generator: DirectedGenerator,
    side_x: Side,
    side_y: Side,
) -> DirectedGenerator:
    """
    Apply R1 rule to a generator by choosing two sides X and Y.
    
    R1 is obtained by choosing two sides X and Y of the generator, such that
    if X = Y then X is not a hybrid vertex, and hanging a new hybrid vertex
    under X and Y.
    
    Parameters
    ----------
    generator : DirectedGenerator
        The level-(k-1) generator to transform.
    side_x : Side
        First side (DirEdgeSide or HybridSide).
    side_y : Side
        Second side (DirEdgeSide or HybridSide).
    
    Returns
    -------
    DirectedGenerator
        Level-k generator resulting from applying R1.
    

    Notes
    -----
    - It is assumed that X and Y are not the same hybrid side.
    - If X = Y (edge side (u, v)): subdivide twice into u -> w1 -> w2 -> v,
      then add edges w1 -> new_hybrid and w2 -> new_hybrid.
    - If X ≠ Y: for edge side (u, v), subdivide into u -> w -> v and add w -> new_hybrid;
      for hybrid side (node h), add edge h -> new_hybrid.
    - The new hybrid node will have in-degree 2 (from the two sides).
    - Level-0 special case: generator has one NodeSide (the single vertex). R1 with
      X = Y = that NodeSide adds two edges from that vertex to a new hybrid (level-1).
    """
    # Level-0 special case: single vertex, attach two edges to new hybrid
    if generator.level == 0:
        if not (side_x == side_y and isinstance(side_x, NodeSide) and not isinstance(side_x, HybridSide)):
            raise PhyloZooValueError(
                "R1 on a level-0 generator requires both sides to be the same NodeSide (IsolatedNodeSide)."
            )
        root = generator.root_node
        new_graph = generator.graph.copy()
        new_hybrid_node = next(new_graph.generate_node_ids(1))
        new_graph.add_node(new_hybrid_node)
        new_graph.add_edge(root, new_hybrid_node)
        new_graph.add_edge(root, new_hybrid_node)
        with no_validation(classes=["DirectedGenerator"]):
            return DirectedGenerator(new_graph)

    # Create a copy of the generator's graph
    new_graph = generator.graph.copy()
    new_hybrid_node = next(new_graph.generate_node_ids(1))
    new_graph.add_node(new_hybrid_node)

    # Case 1: X = Y (same edge side)
    if side_x == side_y and isinstance(side_x, DirEdgeSide):
        # Let uv be X=Y. Subdivide (u,v) twice into (u, w1, w2, v)
        u, v, key = side_x.u, side_x.v, side_x.key
        new_graph.remove_edge(u, v, key=key)
        w1 = next(new_graph.generate_node_ids(1))
        new_graph.add_node(w1)
        w2 = next(new_graph.generate_node_ids(1))
        new_graph.add_node(w2)
        
        # Create path: u -> w1 -> w2 -> v (all directed)
        new_graph.add_edge(u, w1)
        new_graph.add_edge(w1, w2)
        new_graph.add_edge(w2, v)
        
        # Add directed edges: w1 -> new_hybrid and w2 -> new_hybrid
        new_graph.add_edge(w1, new_hybrid_node)
        new_graph.add_edge(w2, new_hybrid_node)
    
    # Case 2: X ≠ Y
    else:
        # Iterate through both sides (they have similar operations)
        for side in [side_x, side_y]:
            if isinstance(side, DirEdgeSide):
                # Subdivide edge (u, v) into u -> w -> v
                u, v, key = side.u, side.v, side.key
                new_graph.remove_edge(u, v, key=key)
                w = next(new_graph.generate_node_ids(1))
                new_graph.add_node(w)
                
                # Create path: u -> w -> v
                new_graph.add_edge(u, w)
                new_graph.add_edge(w, v)
                
                # Add edge: w -> new_hybrid
                new_graph.add_edge(w, new_hybrid_node)
            
            elif isinstance(side, HybridSide):
                # Simply add directed out-edge from hybrid node
                new_graph.add_edge(side.node, new_hybrid_node)
    
    # Create and return the new generator (without validation - will validate later if unique)
    with no_validation(classes=["DirectedGenerator"]):
        new_generator = DirectedGenerator(new_graph)
    return new_generator


def _apply_R2(
    generator: DirectedGenerator,
    side_x: Side,
    side_y: DirEdgeSide,
) -> DirectedGenerator:
    """
    Apply R2 rule to a generator by choosing two sides X and Y.
    
    R2 takes side_x (any side) and side_y (must be an edge side), and adds
    an edge from side_x to side_y.
    
    Parameters
    ----------
    generator : DirectedGenerator
        The level-(k-1) generator to transform.
    side_x : Side
        First side (DirEdgeSide or HybridSide).
    side_y : DirEdgeSide
        Second side (must be an edge side).
    
    Returns
    -------
    DirectedGenerator
        Level-k generator resulting from applying R2.
    
    Notes
    -----
    - side_y must be an edge side (DirEdgeSide).
    - It is assumed that side_x is not reachable from side_y
      (this check is done in the caller before invoking this function).
    - Case 1: Both are edge sides:
      a) If same edge side (u, v): subdivide into u -> w1 -> w2 -> v, then add edge w1 -> w2.
      b) If different edge sides: subdivide both, add edge from subdivision vertex of side_x to subdivision vertex of side_y.
    - Case 2: side_x is hybrid side: add edge from hybrid to subdivision vertex of side_y.
    """
    # Create a copy of the generator's graph
    new_graph = generator.graph.copy()

    # Case 1: Both are edge sides
    if isinstance(side_x, DirEdgeSide) and isinstance(side_y, DirEdgeSide):
        # Case 1a: Same edge side (u, v)
        if side_x == side_y:
            u, v, key = side_x.u, side_x.v, side_x.key
            new_graph.remove_edge(u, v, key=key)
            w1 = next(new_graph.generate_node_ids(1))
            new_graph.add_node(w1)
            w2 = next(new_graph.generate_node_ids(1))
            new_graph.add_node(w2)
            
            # Create path: u -> w1 -> w2 -> v (all directed)
            new_graph.add_edge(u, w1)
            new_graph.add_edge(w1, w2)
            new_graph.add_edge(w2, v)
            
            # Add edge: w1 -> w2
            new_graph.add_edge(w1, w2)
        
        # Case 1b: Different edge sides
        else:
            # Subdivide side_x: (u_x, v_x) -> u_x -> w_x -> v_x
            u_x, v_x, key_x = side_x.u, side_x.v, side_x.key
            new_graph.remove_edge(u_x, v_x, key=key_x)
            w_x = next(new_graph.generate_node_ids(1))
            new_graph.add_node(w_x)
            new_graph.add_edge(u_x, w_x)
            new_graph.add_edge(w_x, v_x)
            # Subdivide side_y: (u_y, v_y) -> u_y -> w_y -> v_y
            u_y, v_y, key_y = side_y.u, side_y.v, side_y.key
            new_graph.remove_edge(u_y, v_y, key=key_y)
            w_y = next(new_graph.generate_node_ids(1))
            new_graph.add_node(w_y)
            new_graph.add_edge(u_y, w_y)
            new_graph.add_edge(w_y, v_y)
            
            # Add edge from subdivision vertex of side_x to subdivision vertex of side_y
            new_graph.add_edge(w_x, w_y)
    
    # Case 2: side_x is hybrid side, side_y is edge side
    elif isinstance(side_x, HybridSide) and isinstance(side_y, DirEdgeSide):
        # Subdivide side_y: (u_y, v_y) -> u_y -> w_y -> v_y
        u_y, v_y, key_y = side_y.u, side_y.v, side_y.key
        new_graph.remove_edge(u_y, v_y, key=key_y)
        w_y = next(new_graph.generate_node_ids(1))
        new_graph.add_node(w_y)
        new_graph.add_edge(u_y, w_y)
        new_graph.add_edge(w_y, v_y)
        
        # Add edge from hybrid node to subdivision vertex of side_y
        new_graph.add_edge(side_x.node, w_y)
    
    # Create and return the new generator (without validation - will validate later if unique)
    with no_validation(classes=["DirectedGenerator"]):
        new_generator = DirectedGenerator(new_graph)
    return new_generator


def _apply_rules(generator: DirectedGenerator) -> Iterator[DirectedGenerator]:
    """
    Apply R1 and R2 rules to all applicable sides of a generator.
    
    This function applies R1 to all pairs of sides and R2 to all hybrid sides,
    yielding all level-k generators that can be constructed from this level-(k-1) generator,
    not up to isomorphism.
    
    Parameters
    ----------
    generator : DirectedGenerator
        The level-(k-1) generator to transform.
    
    Yields
    ------
    DirectedGenerator
        Level-k generators resulting from applying R1 and R2 rules.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
    >>> gen1 = DirectedGenerator(gen_graph)
    >>> level2_generators = list(apply_rules(gen1))
    >>> len(level2_generators) > 0
    True
    """
    sides = generator.sides
    
    # Pre-compute node reachability matrix once (reused for all R2 checks)
    node_reachability = _get_node_reachability_matrix(generator)
    
    def _is_side_reachable_from(side1: Side, side2: Side) -> bool:
        """Check if side1 is reachable from side2 using the node reachability matrix."""
        # Get source node of side1 and target node of side2
        if isinstance(side1, NodeSide):
            source1 = side1.node
        elif isinstance(side1, DirEdgeSide):
            source1 = side1.u
        else:
            return False
        
        if isinstance(side2, NodeSide):
            target2 = side2.node
        elif isinstance(side2, DirEdgeSide):
            target2 = side2.v
        else:
            return False
        
        # Check if source of side1 is reachable from target of side2
        # (i.e., there's a path from target2 to source1)
        return (target2, source1) in node_reachability
    
    # Apply R1 to all pairs of sides (including same side if it's an edge side)
    # Skip pairs where both are hybrid sides and they are the same (X = Y)
    for side_x in sides:
        for side_y in sides:
        # Skip if X = Y and both are hybrid sides
            if side_x == side_y and isinstance(side_x, HybridSide) and isinstance(side_y, HybridSide):
                continue
            
            new_gen = _apply_R1(generator, side_x, side_y)
            yield new_gen
    
    # Apply R2 to all pairs: side_x (any side) and side_y (edge side only),
    # including same sides. Skip if side_x is reachable from side_y.
    edge_sides = generator.edge_sides
    for side_x in sides:
        for side_y in edge_sides:
            # Skip if side_x is reachable from side_y (using helper function)
            # This is mistakenly stated the other way around in the original paper.
            if _is_side_reachable_from(side_x, side_y):
                continue
            
            new_gen = _apply_R2(generator, side_x, side_y)
            yield new_gen



def all_level_k_generators(k: int) -> set[DirectedGenerator]:
    """
    Generate all (strict) level-k generators.
    
    This function constructs all (strict) level-k generators by starting 
    with level-0 generators and iteratively applying R1 and R2 rules.
    
    Parameters
    ----------
    k : int
        The level of generators to generate.
    
    Returns
    -------
    set[DirectedGenerator]
        Set of all level-k generators (up to isomorphism).
    
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
    >>> # Get all level-1 generators
    >>> level1 = all_level_k_generators(1)
    >>> len(level1) == 1
    True
    >>> # Get all level-2 generators
    >>> level2 = all_level_k_generators(2)
    >>> len(level2) == 4
    True
    """
    if k < 0:
        raise PhyloZooValueError("Level must be non-negative")

    if k == 0:
        # Level-0 generators are single nodes
        # Return set with one generator (single node)
        gen_graph = DirectedMultiGraph()
        gen_graph.add_node(0)  # Use 0 as the node ID
        return {DirectedGenerator(gen_graph)}

    # Get all level-(k-1) generators (k>=1: apply R1/R2 to previous level)
    prev_level_generators = all_level_k_generators(k - 1)
    
    # Apply R1 and R2 to each, filtering out isomorphic generators
    # Use invariants to group candidates and reduce isomorphism checks
    result: list[DirectedGenerator] = []
    # Dictionary mapping invariant -> list of generators with that invariant
    invariant_groups: dict[tuple[int, int, tuple[int, ...], tuple[int, ...], tuple[int, ...]], list[DirectedGenerator]] = {}
    
    for prev_gen in prev_level_generators:
        for new_gen in _apply_rules(prev_gen):
            # Compute invariants for fast filtering
            invariant = _get_graph_invariant(new_gen.graph)
            
            # Only check isomorphism against generators with matching invariants
            candidate_group = invariant_groups.get(invariant)
            is_duplicate = False
            
            if candidate_group is not None:
                # Only check against candidates with matching invariants
                for existing_gen in candidate_group:
                    if is_isomorphic(new_gen.graph, existing_gen.graph):
                        is_duplicate = True
                        break
            
            # Only add if not isomorphic to any existing generator
            if not is_duplicate:
                result.append(new_gen)
                # Add to invariant group for future comparisons
                if candidate_group is None:
                    invariant_groups[invariant] = [new_gen]
                else:
                    candidate_group.append(new_gen)

    return set(result)

