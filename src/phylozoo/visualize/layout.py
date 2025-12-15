"""
Layout computation functions for phylogenetic networks.

This module provides optimal layout algorithms for DirectedPhyNetwork and
SemiDirectedPhyNetwork, with crossing minimization and parallel edge handling.
"""

from __future__ import annotations

import random
from typing import TYPE_CHECKING, TypeVar

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork

T = TypeVar('T')


def compute_directed_layout(
    network: 'DirectedPhyNetwork',
    layer_gap: float = 1.5,
    leaf_gap: float = 1.0,
    trials: int = 2000,
    seed: int | None = None,
    orientation: str = 'top-bottom',
    use_branch_lengths: bool = False,
    x_scale: float = 1.5,
    y_scale: float = 1.0,
) -> dict[T, tuple[float, float]]:
    """
    Compute optimal hierarchical layout for DirectedPhyNetwork.
    
    Uses a tree-backbone heuristic with crossing minimization to position
    nodes in a hierarchical layout. Supports top-bottom or left-right
    orientations, with optional branch length scaling.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to layout.
    layer_gap : float, optional
        Spacing between hierarchical layers. By default 1.5.
    leaf_gap : float, optional
        Spacing between leaves within a layer. By default 1.0.
    trials : int, optional
        Number of random child orderings to try for optimization. By default 2000.
    seed : int | None, optional
        Random seed for reproducibility. By default None.
    orientation : str, optional
        Layout direction: 'top-bottom' (root at top) or 'left-right' (root at left).
        By default 'top-bottom'.
    use_branch_lengths : bool, optional
        If True, scale node positions based on branch lengths. By default False.
    x_scale : float, optional
        Scaling factor for x coordinates. By default 1.5.
    y_scale : float, optional
        Scaling factor for y coordinates. By default 1.0.
    
    Returns
    -------
    dict[T, tuple[float, float]]
        Dictionary mapping node IDs to (x, y) positions.
    
    Raises
    ------
    ValueError
        If network is not a DAG or has no root node.
    
    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.layout import compute_directed_layout
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> pos = compute_directed_layout(net)
    >>> len(pos)
    3
    """
    # Get root node
    root = network.root_node
    
    # Create NetworkX DiGraph for layout computation
    G_nx = nx.DiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        G_nx.add_edge(u, v)
    
    if not nx.is_directed_acyclic_graph(G_nx):
        raise ValueError("Network must be a directed acyclic graph (DAG)")
    
    rng = random.Random(seed)
    
    # Step 1: Build tree backbone (spanning tree)
    topo_order = list(nx.topological_sort(G_nx))
    topo_index = {n: i for i, n in enumerate(topo_order)}
    
    tree_edges = []
    for node in topo_order:
        preds = list(G_nx.predecessors(node))
        if preds:
            # Choose parent with earliest topological order
            parent = min(preds, key=lambda p: topo_index[p])
            tree_edges.append((parent, node))
    
    T = nx.DiGraph(tree_edges)
    
    # Step 2: Recursive layout helper (x-coordinate assignment)
    def layout_tree(node: T, depth: int, x_offset: list[float]) -> tuple[dict[T, tuple[float, float]], float]:
        """Recursively layout tree, returning positions and next x offset."""
        children = list(T.successors(node))
        if not children:
            # Leaf node
            return {node: (x_offset[0], depth)}, x_offset[0] + leaf_gap
        
        pos = {}
        child_xs = []
        for ch in children:
            child_pos, x_next = layout_tree(ch, depth + 1, x_offset)
            pos.update(child_pos)
            child_xs.append(pos[ch][0])
            x_offset[0] = x_next
        
        # Position parent at mean of children's x-coordinates
        mean_x = np.mean(child_xs) if child_xs else x_offset[0]
        pos[node] = (mean_x, depth)
        return pos, x_offset[0]
    
    # Step 3: Crossing count heuristic
    def count_crossings(
        pos: dict[T, tuple[float, float]],
        edges: list[tuple[T, T]]
    ) -> int:
        """Count number of edge crossings."""
        xs = [(pos[u][0], pos[v][0]) for u, v in edges if u in pos and v in pos]
        count = 0
        for i in range(len(xs)):
            x1u, x1v = xs[i]
            for j in range(i + 1, len(xs)):
                x2u, x2v = xs[j]
                # Check if edges cross (x-coordinates are interleaved)
                if (x1u - x2u) * (x1v - x2v) < 0:
                    count += 1
        return count
    
    # Step 4: Optimization loop (minimize crossings)
    best_pos: dict[T, tuple[float, float]] | None = None
    best_score = float('inf')
    node_children = {n: list(T.successors(n)) for n in T.nodes}
    
    for _ in range(trials):
        # Randomly shuffle child orderings
        for n in node_children:
            rng.shuffle(node_children[n])
        
        # Rebuild tree with new ordering
        T_tmp = nx.DiGraph()
        for u in node_children:
            for v in node_children[u]:
                T_tmp.add_edge(u, v)
        
        # Compute layout
        pos, _ = layout_tree(root, 0, [0.0])
        
        # Count crossings (tree edges + extra edges)
        extra_edges = [e for e in G_nx.edges if e not in tree_edges]
        score = count_crossings(pos, list(tree_edges) + extra_edges)
        
        if score < best_score:
            best_pos = pos.copy()
            best_score = score
    
    if best_pos is None:
        raise ValueError("Failed to compute layout")
    
    pos = best_pos
    
    # Step 5: Assign depth coordinates
    depths: dict[T, int] = {}
    for node in topo_order:
        preds = list(G_nx.predecessors(node))
        if preds:
            depths[node] = max(depths[p] for p in preds) + 1
        else:
            depths[node] = 0
    max_depth = max(depths.values()) if depths else 0
    
    # Step 6: Apply branch lengths if requested
    if use_branch_lengths:
        pos = _apply_branch_lengths(network, pos, depths, max_depth, orientation)
    
    # Step 7: Final coordinate mapping
    if orientation.lower() == 'top-bottom':
        # Root at top, children below
        final_pos = {}
        for n in pos:
            x, _ = pos[n]
            y = (max_depth - depths[n]) * layer_gap
            final_pos[n] = (x * x_scale, y * y_scale)
        return final_pos
    elif orientation.lower() == 'left-right':
        # Root on left, children to the right
        final_pos = {}
        for n in pos:
            y, _ = pos[n]  # Reuse previous 'x' as vertical coordinate
            x = depths[n] * layer_gap
            final_pos[n] = (x * x_scale, y * y_scale)
        return final_pos
    else:
        raise ValueError(f"orientation must be 'top-bottom' or 'left-right', got '{orientation}'")


def compute_semidirected_layout(
    network: 'MixedPhyNetwork | SemiDirectedPhyNetwork',
    iterations: int = 50,
    k: float = 2.0,
    seed: int | None = None,
    leaf_repulsion: float = 1.5,
) -> dict[T, tuple[float, float]]:
    """
    Compute force-directed layout for SemiDirectedPhyNetwork with leaf repulsion.
    
    Uses a force-directed algorithm with extra repulsion for leaf nodes to push
    them outward, creating an unrooted tree-like appearance.
    
    Parameters
    ----------
    network : MixedPhyNetwork | SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to layout.
    iterations : int, optional
        Number of iterations for force-directed algorithm. By default 50.
    k : float, optional
        Optimal distance between nodes. By default 2.0.
    seed : int | None, optional
        Random seed for reproducibility. By default None.
    leaf_repulsion : float, optional
        Multiplier for leaf repulsion force. Higher values push leaves further out.
        By default 1.5.
    
    Returns
    -------
    dict[T, tuple[float, float]]
        Dictionary mapping node IDs to (x, y) positions.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.visualize.layout import compute_semidirected_layout
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> pos = compute_semidirected_layout(net)
    >>> len(pos)
    3
    """
    # Get leaves
    leaves = network.leaves
    
    # Create NetworkX MultiGraph for layout
    G_nx = nx.MultiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    
    # Add all edges (directed and undirected)
    for edge_tuple in network._graph.directed_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            G_nx.add_edge(u, v, key=key)
    for edge_tuple in network._graph.undirected_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            G_nx.add_edge(u, v, key=key)
    
    # Start with spring layout
    pos = nx.spring_layout(G_nx, k=k, iterations=iterations, seed=seed)
    
    # Apply leaf repulsion
    pos = _apply_leaf_repulsion(pos, leaves, leaf_repulsion, iterations)
    
    return pos


def _apply_branch_lengths(
    network: 'DirectedPhyNetwork',
    pos: dict[T, tuple[float, float]],
    depths: dict[T, int],
    max_depth: int,
    orientation: str,
) -> dict[T, tuple[float, float]]:
    """
    Apply branch length scaling to node positions.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The network with branch lengths.
    pos : dict[T, tuple[float, float]]
        Current node positions.
    depths : dict[T, int]
        Depth of each node.
    max_depth : int
        Maximum depth in the network.
    orientation : str
        Layout orientation ('top-bottom' or 'left-right').
    
    Returns
    -------
    dict[T, tuple[float, float]]
        Updated positions with branch length scaling.
    """
    # Collect all branch lengths
    branch_lengths: list[float] = []
    for u, v, key in network._graph.edges(keys=True):
        bl = network.get_branch_length(u, v, key)
        if bl is not None and bl > 0:
            branch_lengths.append(bl)
    
    if not branch_lengths:
        return pos
    
    # Normalize branch lengths to reasonable range (0.5 to 2.0)
    min_bl = min(branch_lengths)
    max_bl = max(branch_lengths)
    if max_bl == min_bl:
        scale = 1.0
    else:
        # Map to [0.5, 2.0] range
        scale_range = 1.5  # 2.0 - 0.5
        scale = scale_range / (max_bl - min_bl)
    
    # Update positions based on branch lengths
    new_pos = pos.copy()
    root = network.root_node
    
    # BFS from root, adjusting positions
    queue = [(root, pos[root])]
    visited = {root}
    
    while queue:
        node, (x, y) = queue.pop(0)
        
        for child in network.children(node):
            if child in visited:
                continue
            visited.add(child)
            
            # Get branch length (use first edge if multiple parallel edges)
            bl = network.get_branch_length(node, child)
            if bl is None or bl <= 0:
                bl = 1.0  # Default length
            
            # Normalize branch length
            normalized_bl = 0.5 + (bl - min_bl) * scale if max_bl != min_bl else 1.0
            
            # Adjust child position
            child_x, child_y = pos[child]
            if orientation.lower() == 'top-bottom':
                # Adjust y-coordinate (vertical distance)
                new_y = y - normalized_bl
                new_pos[child] = (child_x, new_y)
            else:  # left-right
                # Adjust x-coordinate (horizontal distance)
                new_x = x + normalized_bl
                new_pos[child] = (new_x, child_y)
            
            queue.append((child, new_pos[child]))
    
    return new_pos


def _apply_leaf_repulsion(
    pos: dict[T, tuple[float, float]],
    leaves: set[T],
    repulsion_strength: float,
    iterations: int,
) -> dict[T, tuple[float, float]]:
    """
    Apply repulsion force to push leaf nodes outward.
    
    Parameters
    ----------
    pos : dict[T, tuple[float, float]]
        Current node positions.
    leaves : set[T]
        Set of leaf nodes.
    repulsion_strength : float
        Strength of repulsion force.
    iterations : int
        Number of iterations to apply repulsion.
    
    Returns
    -------
    dict[T, tuple[float, float]]
        Updated positions with leaves pushed outward.
    """
    if not pos or not leaves:
        return pos
    
    # Compute center of mass
    all_x = [x for x, y in pos.values()]
    all_y = [y for x, y in pos.values()]
    center_x = np.mean(all_x) if all_x else 0.0
    center_y = np.mean(all_y) if all_y else 0.0
    
    new_pos = pos.copy()
    
    # Apply repulsion iteratively
    for _ in range(iterations):
        for node in leaves:
            if node not in new_pos:
                continue
            
            x, y = new_pos[node]
            
            # Vector from center to node
            dx = x - center_x
            dy = y - center_y
            dist = np.sqrt(dx**2 + dy**2) if (dx != 0 or dy != 0) else 1.0
            
            if dist > 0:
                # Normalize and apply repulsion
                dx_norm = dx / dist
                dy_norm = dy / dist
                
                # Repulsion force (push away from center)
                force = repulsion_strength * 0.01  # Small step size
                new_pos[node] = (
                    x + dx_norm * force,
                    y + dy_norm * force
                )
    
    return new_pos


def _count_crossings(
    pos: dict[T, tuple[float, float]],
    edges: list[tuple[T, T]],
) -> int:
    """
    Count the number of edge crossings in a layout.
    
    Parameters
    ----------
    pos : dict[T, tuple[float, float]]
        Node positions.
    edges : list[tuple[T, T]]
        List of edges as (u, v) tuples.
    
    Returns
    -------
    int
        Number of edge crossings.
    """
    crossings = 0
    for i in range(len(edges)):
        u1, v1 = edges[i]
        if u1 not in pos or v1 not in pos:
            continue
        x1u, y1u = pos[u1]
        x1v, y1v = pos[v1]
        
        for j in range(i + 1, len(edges)):
            u2, v2 = edges[j]
            if u2 not in pos or v2 not in pos:
                continue
            x2u, y2u = pos[u2]
            x2v, y2v = pos[v2]
            
            # Check if edges cross (simplified 2D line intersection)
            if _edges_cross((x1u, y1u), (x1v, y1v), (x2u, y2u), (x2v, y2v)):
                crossings += 1
    
    return crossings


def _edges_cross(
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
    p4: tuple[float, float],
) -> bool:
    """
    Check if two line segments cross.
    
    Parameters
    ----------
    p1, p2 : tuple[float, float]
        Endpoints of first edge.
    p3, p4 : tuple[float, float]
        Endpoints of second edge.
    
    Returns
    -------
    bool
        True if edges cross, False otherwise.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    
    # Check if segments are parallel
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(denom) < 1e-10:
        return False
    
    # Compute intersection point
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom
    
    # Check if intersection is within both segments
    return 0 <= t <= 1 and 0 <= u <= 1

