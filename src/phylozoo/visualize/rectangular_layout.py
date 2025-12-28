"""
Rectangular layout computation for directed phylogenetic networks.

This module provides a Dendroscope-style rectangular layout algorithm that
positions nodes in strict horizontal layers with crossing minimization.
Based on the approach described in Huson (2025) "Sketch, capture and layout phylogenies".
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


def compute_rectangular_layout(
    network: 'DirectedPhyNetwork',
    layer_gap: float = 1.5,
    node_gap: float = 1.0,
    iterations: int = 50,
    seed: int | None = None,
    orientation: str = 'top-bottom',
    use_branch_lengths: bool = False,
) -> dict[T, tuple[float, float]]:
    """
    Compute rectangular layout for DirectedPhyNetwork (Dendroscope-style).
    
    This layout algorithm positions nodes in strict horizontal layers based on
    their distance from the root, and uses a barycenter heuristic to minimize
    edge crossings within layers. This produces a "combining view" layout similar
    to Dendroscope and ape.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to layout.
    layer_gap : float, optional
        Vertical spacing between layers. By default 1.5.
    node_gap : float, optional
        Horizontal spacing between nodes within a layer. By default 1.0.
    iterations : int, optional
        Number of iterations for crossing minimization. By default 50.
    seed : int | None, optional
        Random seed for reproducibility. By default None.
    orientation : str, optional
        Layout direction: 'top-bottom' (root at top) or 'left-right' (root at left).
        By default 'top-bottom'.
    use_branch_lengths : bool, optional
        If True, scale layer positions based on branch lengths. By default False.
    
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
    >>> from phylozoo.visualize.rectangular_layout import compute_rectangular_layout
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> pos = compute_rectangular_layout(net)
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
    
    # Step 1: Assign nodes to layers based on distance from root
    layers = _assign_layers(G_nx, root)
    max_layer = max(layers.values()) if layers else 0
    
    # Step 2: Minimize crossings using barycenter heuristic
    layer_orderings = _minimize_crossings_sugiyama(
        G_nx, layers, max_layer, iterations, seed
    )
    
    # Step 3: Assign x-coordinates within each layer
    pos = {}
    for layer_num in range(max_layer + 1):
        nodes_in_layer = [n for n, l in layers.items() if l == layer_num]
        if not nodes_in_layer:
            continue
        
        # Sort nodes by their ordering
        nodes_in_layer.sort(key=lambda n: layer_orderings.get(n, 0))
        
        # Assign x-coordinates with equal spacing
        num_nodes = len(nodes_in_layer)
        if num_nodes == 1:
            x_positions = [0.0]
        else:
            total_width = (num_nodes - 1) * node_gap
            x_positions = np.linspace(-total_width / 2, total_width / 2, num_nodes)
        
        for i, node in enumerate(nodes_in_layer):
            x = x_positions[i]
            y = layer_num * layer_gap
            pos[node] = (x, y)
    
    # Step 4: Apply branch lengths if requested
    if use_branch_lengths:
        pos = _apply_branch_lengths_rectangular(network, pos, layers, max_layer, orientation)
    
    # Step 5: Final coordinate mapping based on orientation
    if orientation.lower() == 'top-bottom':
        # Root at top, children below (y increases downward)
        final_pos = {}
        for n in pos:
            x, y = pos[n]
            final_pos[n] = (x, max_layer * layer_gap - y)
        return final_pos
    elif orientation.lower() == 'left-right':
        # Root on left, children to the right (x increases rightward)
        final_pos = {}
        for n in pos:
            y, x = pos[n]  # Swap coordinates
            final_pos[n] = (x, y)
        return final_pos
    else:
        raise ValueError(f"orientation must be 'top-bottom' or 'left-right', got '{orientation}'")


def _assign_layers(G_nx: nx.DiGraph, root: T) -> dict[T, int]:
    """
    Assign nodes to layers based on longest path from root.
    
    Uses topological sort to ensure all predecessors are processed before
    assigning a layer to a node. The layer is the maximum distance from root.
    After initial assignment, all leaves are moved to the bottom layer to create
    a more Dendroscope-like rectangular layout.
    
    Parameters
    ----------
    G_nx : nx.DiGraph
        NetworkX directed graph representation.
    root : T
        Root node of the network.
    
    Returns
    -------
    dict[T, int]
        Dictionary mapping node to layer number (0 = root layer).
    """
    layers: dict[T, int] = {root: 0}
    
    # Use topological sort to process nodes in correct order
    topo_order = list(nx.topological_sort(G_nx))
    
    # First pass: assign layers based on longest path
    for node in topo_order:
        if node == root:
            continue
        
        # Get all predecessors
        predecessors = list(G_nx.predecessors(node))
        if not predecessors:
            # Node with no predecessors (shouldn't happen in DAG with root)
            layers[node] = 0
        else:
            # Layer is max of all predecessor layers + 1
            max_pred_layer = max(layers.get(p, 0) for p in predecessors)
            layers[node] = max_pred_layer + 1
    
    # Second pass: identify leaves and find max layer
    leaves = [n for n in G_nx.nodes if G_nx.out_degree(n) == 0]
    max_layer = max(layers.values()) if layers else 0
    
    # Move all leaves to the bottom layer (max_layer)
    # This creates a more Dendroscope-like layout where all leaves are at the same level
    for leaf in leaves:
        layers[leaf] = max_layer
    
    # Recompute layers for nodes that are now "above" leaves
    # We need to ensure no node is at a layer >= max_layer unless it's a leaf
    changed = True
    while changed:
        changed = False
        for node in topo_order:
            if node == root or node in leaves:
                continue
            
            # Get all predecessors
            predecessors = list(G_nx.predecessors(node))
            if predecessors:
                max_pred_layer = max(layers.get(p, 0) for p in predecessors)
                new_layer = max_pred_layer + 1
                # Ensure node is not at or below leaf layer (unless it's a leaf)
                if new_layer >= max_layer:
                    new_layer = max_layer - 1
                if layers.get(node, 0) != new_layer:
                    layers[node] = new_layer
                    changed = True
    
    return layers


def _minimize_crossings_sugiyama(
    G_nx: nx.DiGraph,
    layers: dict[T, int],
    max_layer: int,
    iterations: int,
    seed: int | None,
) -> dict[T, float]:
    """
    Minimize edge crossings using barycenter heuristic (Sugiyama-style).
    
    Parameters
    ----------
    G_nx : nx.DiGraph
        NetworkX directed graph representation.
    layers : dict[T, int]
        Layer assignment for each node.
    max_layer : int
        Maximum layer number.
    iterations : int
        Number of iterations for refinement.
    seed : int | None
        Random seed.
    
    Returns
    -------
    dict[T, float]
        Dictionary mapping node to its ordering value (for sorting within layer).
    """
    import random
    rng = random.Random(seed)
    
    # Initialize orderings randomly
    orderings: dict[T, float] = {}
    for node in G_nx.nodes:
        orderings[node] = rng.random()
    
    # Iteratively refine using barycenter heuristic
    for _ in range(iterations):
        new_orderings: dict[T, float] = {}
        
        # Process layers from top to bottom
        for layer_num in range(max_layer + 1):
            nodes_in_layer = [n for n, l in layers.items() if l == layer_num]
            
            for node in nodes_in_layer:
                # Compute barycenter from parents (nodes in previous layer)
                parents = list(G_nx.predecessors(node))
                if parents:
                    # Average x-position of parents
                    parent_positions = [orderings.get(p, 0.0) for p in parents]
                    new_orderings[node] = np.mean(parent_positions) if parent_positions else orderings.get(node, 0.0)
                else:
                    # Root or node with no parents: keep current or random
                    new_orderings[node] = orderings.get(node, rng.random())
        
        # Update orderings
        orderings.update(new_orderings)
    
    return orderings


def _apply_branch_lengths_rectangular(
    network: 'DirectedPhyNetwork',
    pos: dict[T, tuple[float, float]],
    layers: dict[T, int],
    max_layer: int,
    orientation: str,
) -> dict[T, tuple[float, float]]:
    """
    Apply branch length scaling to rectangular layout.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The network with branch lengths.
    pos : dict[T, tuple[float, float]]
        Current node positions.
    layers : dict[T, int]
        Layer assignment for each node.
    max_layer : int
        Maximum layer number.
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
    
    # Normalize branch lengths to reasonable range (0.5 to 2.0 times layer_gap)
    min_bl = min(branch_lengths)
    max_bl = max(branch_lengths)
    if max_bl == min_bl:
        scale = 1.0
    else:
        scale_range = 1.5  # 2.0 - 0.5
        scale = scale_range / (max_bl - min_bl)
    
    # Update positions based on branch lengths
    new_pos = pos.copy()
    root = network.root_node
    
    # BFS from root, adjusting layer positions
    queue = [(root, pos[root])]
    visited = {root}
    
    while queue:
        node, (x, y) = queue.pop(0)
        
        for child in network.children(node):
            if child in visited:
                continue
            visited.add(child)
            
            # Get branch length
            bl = network.get_branch_length(node, child)
            if bl is None or bl <= 0:
                bl = 1.0  # Default length
            
            # Normalize branch length
            normalized_bl = 0.5 + (bl - min_bl) * scale if max_bl != min_bl else 1.0
            
            # Adjust child position
            child_x, child_y = pos[child]
            if orientation.lower() == 'top-bottom':
                # Adjust y-coordinate (vertical distance)
                new_y = y + normalized_bl
                new_pos[child] = (child_x, new_y)
            else:  # left-right
                # Adjust x-coordinate (horizontal distance)
                new_x = x + normalized_bl
                new_pos[child] = (new_x, child_y)
            
            queue.append((child, new_pos[child]))
    
    return new_pos

