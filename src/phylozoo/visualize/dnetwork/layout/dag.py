"""
DAG layout algorithm for DirectedPhyNetwork.

This module implements a tree-backbone heuristic layout with crossing minimization.
Based on a tree-backbone approach that optimizes child ordering to minimize edge crossings.
"""

from __future__ import annotations

import random
from typing import TYPE_CHECKING, TypeVar

import networkx as nx
import numpy as np

from .base import DAGLayout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


def compute_dag_layout(
    network: 'DirectedPhyNetwork',
    layer_gap: float = 1.5,
    leaf_gap: float = 1.0,
    trials: int = 2000,
    seed: int | None = None,
    direction: str = 'TD',  # "TD" = top-down, "LR" = left-right
    x_scale: float = 1.5,
    y_scale: float = 1.0,
) -> DAGLayout:
    """
    Layout a phylogenetic network (DAG) using a tree-backbone heuristic.

    This function implements a tree-backbone layout algorithm:
    1. Builds a tree backbone (one parent per node)
    2. Recursively layouts the tree (x coordinates based on children)
    3. Optimizes child ordering to minimize edge crossings
    4. Assigns depth coordinates
    5. Maps to final coordinates based on direction

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    layer_gap : float, optional
        Spacing between hierarchical layers. By default 1.5.
    leaf_gap : float, optional
        Spacing between leaves within a layer. By default 1.0.
    trials : int, optional
        Number of random child orderings to try for optimization.
        By default 2000.
    seed : int | None, optional
        Random seed for reproducibility. By default None.
    direction : str, optional
        Layout direction: 'TD' (top-down) or 'LR' (left-right).
        By default 'TD'.
    x_scale : float, optional
        Scaling factor applied to x coordinates. By default 1.5.
    y_scale : float, optional
        Scaling factor applied to y coordinates. By default 1.0.

    Returns
    -------
    DAGLayout
        The computed layout.

    Raises
    ------
    ValueError
        If network is empty, not a DAG, or direction is invalid.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.dnetwork.layout import compute_dag_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_dag_layout(net)
    >>> len(layout.positions)
    3
    """
    if network.number_of_nodes() == 0:
        raise ValueError("Cannot compute layout for empty network")

    if direction.upper() not in ('TD', 'LR'):
        raise ValueError(f"direction must be 'TD' or 'LR', got '{direction}'")

    rng = random.Random(seed)

    # Convert network to NetworkX DiGraph for layout computation
    G = nx.DiGraph()
    for node in network._graph.nodes:
        G.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        G.add_edge(u, v)

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("Network must be a DAG")

    # --- Step 1: Build a tree backbone ---
    roots = [n for n in G.nodes if G.in_degree(n) == 0]
    if not roots:
        raise ValueError("No root found (no node with in-degree 0).")
    root = roots[0]

    topo_order = list(nx.topological_sort(G))
    topo_index = {n: i for i, n in enumerate(topo_order)}

    tree_edges: list[tuple[T, T]] = []
    for node in topo_order:
        preds = list(G.predecessors(node))
        if preds:
            parent = min(preds, key=lambda p: topo_index[p])
            tree_edges.append((parent, node))

    T = nx.DiGraph(tree_edges)

    # --- Step 2: Recursive layout helper (x only) ---
    def layout_tree(node: T, depth: int, x_offset: list[float]) -> tuple[dict[T, tuple[float, float]], float]:
        """Recursively layout tree, returning positions and next x offset."""
        children = list(T.successors(node))
        if not children:
            pos = {node: (x_offset[0], depth)}
            return pos, x_offset[0] + leaf_gap

        pos: dict[T, tuple[float, float]] = {}
        child_xs: list[float] = []
        for ch in children:
            child_pos, x_next = layout_tree(ch, depth + 1, x_offset)
            pos.update(child_pos)
            child_xs.append(pos[ch][0])
            x_offset[0] = x_next

        mean_x = float(np.mean(child_xs))
        pos[node] = (mean_x, depth)
        return pos, x_offset[0]

    # --- Step 3: Crossing count heuristic ---
    def count_crossings(
        pos: dict[T, tuple[float, float]], edges: list[tuple[T, T]]
    ) -> int:
        """Count edge crossings."""
        xs = [
            (pos[u][0], pos[v][0])
            for u, v in edges
            if u in pos and v in pos
        ]
        count = 0
        for i in range(len(xs)):
            x1u, x1v = xs[i]
            for j in range(i + 1, len(xs)):
                x2u, x2v = xs[j]
                if (x1u - x2u) * (x1v - x2v) < 0:
                    count += 1
        return count

    # --- Step 4: Optimization loop (for x-ordering) ---
    best_pos: dict[T, tuple[float, float]] | None = None
    best_score = float('inf')
    node_children: dict[T, list[T]] = {n: list(T.successors(n)) for n in T.nodes}

    for _ in range(trials):
        # Shuffle children for this trial
        for n in node_children:
            rng.shuffle(node_children[n])

        # Build temporary tree with shuffled children
        T_tmp = nx.DiGraph()
        for u in node_children:
            for v in node_children[u]:
                T_tmp.add_edge(u, v)

        # Compute layout with this ordering
        pos, _ = layout_tree(root, 0, [0.0])
        
        # Count crossings including non-tree edges
        all_edges = [(u, v) for u, v in G.edges]
        score = count_crossings(pos, all_edges)

        if score < best_score:
            best_pos, best_score = pos, score

    if best_pos is None:
        raise ValueError("Failed to compute layout")

    pos = best_pos

    # --- Step 5: Assign depth coordinates ---
    depths: dict[T, int] = {}
    for node in topo_order:
        preds = list(G.predecessors(node))
        if preds:
            depths[node] = max(depths[p] for p in preds) + 1
        else:
            depths[node] = 0
    max_depth = max(depths.values()) if depths else 0

    # --- Step 6: Final coordinate mapping ---
    final_positions: dict[T, tuple[float, float]] = {}
    if direction.upper() == 'TD':
        # Root at top, children below
        for n in pos:
            x, _ = pos[n]
            y = (max_depth - depths[n]) * layer_gap
            final_positions[n] = (x * x_scale, y * y_scale)
    elif direction.upper() == 'LR':
        # Root on left, children to the right
        for n in pos:
            y, _ = pos[n]  # reuse previous 'x' as vertical coordinate
            x = depths[n] * layer_gap
            final_positions[n] = (x * x_scale, y * y_scale)

    # --- Step 7: Compute edge routes ---
    from ..rendering.routes import compute_backbone_routes, compute_hybrid_routes

    # Identify tree edges and hybrid edges
    tree_edge_set: set[tuple[T, T, int]] = set()
    hybrid_edge_set: set[tuple[T, T, int]] = set()

    for u, v, key in network._graph.edges(keys=True):
        if (u, v) in tree_edges:
            tree_edge_set.add((u, v, key))
        else:
            hybrid_edge_set.add((u, v, key))

    # Compute routes
    backbone_routes = compute_backbone_routes(
        network, final_positions, tree_edge_set
    )
    hybrid_routes = compute_hybrid_routes(
        network,
        final_positions,
        hybrid_edge_set,
    )

    all_routes = {**backbone_routes, **hybrid_routes}

    # --- Step 8: Create layout object ---
    return DAGLayout(
        network=network,
        positions=final_positions,
        edge_routes=all_routes,
        backbone_edges=tree_edge_set,
        reticulate_edges=hybrid_edge_set,
        algorithm='dag',
        parameters={
            'layer_gap': layer_gap,
            'leaf_gap': leaf_gap,
            'trials': trials,
            'seed': seed,
            'direction': direction,
            'x_scale': x_scale,
            'y_scale': y_scale,
        },
    )

