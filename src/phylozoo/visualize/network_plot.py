"""
Network plotting functions.

This module provides functions for plotting and visualizing phylogenetic networks.
"""

from __future__ import annotations

from typing import Any, TypeVar

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np

from .layout import compute_directed_layout, compute_semidirected_layout

# Try to import pyqtgraph (optional dependency)
try:
    import pyqtgraph as pg
    HAS_PYQTGAPH = True
except ImportError:
    HAS_PYQTGAPH = False

T = TypeVar('T')


def plot_network(
    network: 'DirectedPhyNetwork | MixedPhyNetwork | SemiDirectedPhyNetwork',
    ax: plt.Axes | None = None,
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 500,
    leaf_size: int = 600,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    label_offset: float = 0.1,
    orientation: str = 'top-bottom',
    **kwargs
) -> plt.Axes:
    """
    Plot a phylogenetic network (DirectedPhyNetwork, MixedPhyNetwork, or SemiDirectedPhyNetwork).
    
    For DirectedPhyNetwork: Uses hierarchical layout (top-to-bottom or left-to-right) with leaves
    on one end. For SemiDirectedPhyNetwork: Uses planar layout with labels on outside.
    
    Parameters
    ----------
    network : DirectedPhyNetwork | MixedPhyNetwork | SemiDirectedPhyNetwork
        The phylogenetic network to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, computes layout automatically.
    node_color : str, optional
        Color for internal nodes. Default is 'lightblue'.
    leaf_color : str, optional
        Color for leaf nodes. Default is 'lightgreen'.
    hybrid_color : str, optional
        Color for hybrid nodes. Default is 'salmon'.
    node_size : int, optional
        Size of internal nodes. Default is 500.
    leaf_size : int, optional
        Size of leaf nodes. Default is 600.
    edge_color : str, optional
        Color for regular edges. Default is 'gray'.
    hybrid_edge_color : str, optional
        Color for hybrid edges. Default is 'red'.
    edge_width : float, optional
        Width of edges. Default is 2.0.
    arrow_size : int, optional
        Size of arrows. Default is 20.
    with_labels : bool, optional
        Whether to show node labels. Default is True.
    label_offset : float, optional
        Offset for labels from nodes. Default is 0.1.
    orientation : str, optional
        Layout orientation: 'top-bottom' or 'left-right'. Default is 'top-bottom'.
    **kwargs
        Additional keyword arguments passed to matplotlib functions.
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    
    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.network_plot import plot_network
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> plot_network(net)
    <Axes: >
    """
    # Determine network type and call appropriate function
    network_type = type(network).__name__
    
    if network_type == 'DirectedPhyNetwork':
        return plot_directed_network(
            network, ax=ax, pos=pos, node_color=node_color, leaf_color=leaf_color,
            hybrid_color=hybrid_color, node_size=node_size, leaf_size=leaf_size,
            edge_color=edge_color, hybrid_edge_color=hybrid_edge_color,
            edge_width=edge_width, arrow_size=arrow_size, with_labels=with_labels,
            label_offset=label_offset, orientation=orientation, **kwargs
        )
    elif network_type in ('MixedPhyNetwork', 'SemiDirectedPhyNetwork'):
        return plot_semi_directed_network(
            network, ax=ax, pos=pos, node_color=node_color, leaf_color=leaf_color,
            hybrid_color=hybrid_color, node_size=node_size, leaf_size=leaf_size,
            edge_color=edge_color, hybrid_edge_color=hybrid_edge_color,
            edge_width=edge_width, arrow_size=arrow_size, with_labels=with_labels,
            label_offset=label_offset, **kwargs
        )
    else:
        raise ValueError(f"Unsupported network type: {network_type}")


def plot_directed_network(
    network: 'DirectedPhyNetwork',
    ax: plt.Axes | None = None,
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 500,
    leaf_size: int = 600,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    label_offset: float = 0.1,
    orientation: str = 'top-bottom',
    **kwargs
) -> plt.Axes:
    """
    Plot a DirectedPhyNetwork with hierarchical layout.
    
    Uses hierarchical positioning with root at top (or left) and leaves at bottom (or right).
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, computes hierarchical layout.
    orientation : str, optional
        'top-bottom' (default) or 'left-right'.
    **kwargs
        Additional keyword arguments (see plot_network for full list).
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    
    # Get network properties
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Create NetworkX DiGraph for layout
    G_nx = nx.DiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        G_nx.add_edge(u, v)
    
    # Compute hierarchical layout if not provided
    if pos is None:
        if orientation == 'top-bottom':
            pos = _hierarchical_layout_top_bottom(network, G_nx)
        else:  # left-right
            pos = _hierarchical_layout_left_right(network, G_nx)
    
    # Group edges by (u, v) for parallel edge handling
    edge_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for u, v, key in network._graph.edges(keys=True):
        edge_key = (u, v)
        if edge_key not in edge_groups:
            edge_groups[edge_key] = []
        edge_groups[edge_key].append((u, v, key))
    
    # Draw edges with curvature for parallel edges
    for (u, v), edges in edge_groups.items():
        is_hybrid_edge = v in hybrid_nodes
        edge_col = hybrid_edge_color if is_hybrid_edge else edge_color
        
        num_parallel = len(edges)
        
        if num_parallel == 1:
            # Single edge - draw straight
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=edge_col, width=edge_width,
                arrows=True, arrowsize=arrow_size,
                arrowstyle='->', alpha=0.7, **kwargs
            )
        else:
            # Multiple parallel edges - draw with curvature
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                # Create curved path
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=edge_col, linewidth=edge_width,
                    arrowstyle='->', mutation_scale=arrow_size,
                    alpha=0.7, zorder=1
                )
                ax.add_patch(arrow)
    
    # Draw nodes with different colors for leaves and hybrids
    node_colors_map = {}
    node_sizes_map = {}
    for node in network._graph.nodes:
        if node in leaves:
            node_colors_map[node] = leaf_color
            node_sizes_map[node] = leaf_size
        elif node in hybrid_nodes:
            node_colors_map[node] = hybrid_color
            node_sizes_map[node] = node_size
        else:
            node_colors_map[node] = node_color
            node_sizes_map[node] = node_size
    
    nx.draw_networkx_nodes(
        G_nx, pos, ax=ax, node_color=[node_colors_map[n] for n in G_nx.nodes],
        node_size=[node_sizes_map[n] for n in G_nx.nodes], alpha=0.9, **kwargs
    )
    
    # Draw labels
    if with_labels:
        labels = {}
        for node in network._graph.nodes:
            label = network.get_label(node)
            if label:
                labels[node] = label
            else:
                labels[node] = str(node)
        
        # Position labels offset from nodes
        label_pos = {}
        if orientation == 'top-bottom':
            for node, (x, y) in pos.items():
                if node in leaves:
                    label_pos[node] = (x, y - label_offset)
                else:
                    label_pos[node] = (x, y + label_offset)
        else:  # left-right
            for node, (x, y) in pos.items():
                if node in leaves:
                    label_pos[node] = (x + label_offset, y)
                else:
                    label_pos[node] = (x - label_offset, y)
        
        nx.draw_networkx_labels(G_nx, label_pos, labels=labels, ax=ax, font_size=10, font_weight='bold')
    
    ax.set_title('Directed Phylogenetic Network', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return ax


def plot_semi_directed_network(
    network: 'MixedPhyNetwork | SemiDirectedPhyNetwork',
    ax: plt.Axes | None = None,
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 500,
    leaf_size: int = 600,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    label_offset: float = 0.15,
    **kwargs
) -> plt.Axes:
    """
    Plot a MixedPhyNetwork or SemiDirectedPhyNetwork with planar layout.
    
    Uses planar/spring layout with labels positioned on the outside.
    
    Parameters
    ----------
    network : MixedPhyNetwork | SemiDirectedPhyNetwork
        The mixed/semi-directed phylogenetic network to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, computes planar layout.
    **kwargs
        Additional keyword arguments (see plot_network for full list).
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    
    # Get network properties
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Create NetworkX MultiGraph for layout
    G_nx = nx.MultiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    # Add all edges for layout
    for edge_tuple in network._graph.directed_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            G_nx.add_edge(u, v, key=key)
    for edge_tuple in network._graph.undirected_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            G_nx.add_edge(u, v, key=key)
    
    # Compute layout if not provided
    if pos is None:
        # Try planar layout first, fall back to spring layout
        try:
            pos = nx.planar_layout(G_nx)
        except nx.NetworkXException:
            pos = nx.spring_layout(G_nx, seed=42, k=2.0, iterations=50)
    
    # Draw undirected edges with curvature for parallel edges
    undirected_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for edge_tuple in network._graph.undirected_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            edge_key = (min(u, v), max(u, v))
            if edge_key not in undirected_groups:
                undirected_groups[edge_key] = []
            undirected_groups[edge_key].append((u, v, key))
    
    for (u, v), edges in undirected_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=edge_color, width=edge_width,
                style='dashed', alpha=0.7, **kwargs
            )
        else:
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=edge_color, linewidth=edge_width,
                    linestyle='dashed', alpha=0.7, zorder=1,
                    arrowstyle='-'
                )
                ax.add_patch(arrow)
    
    # Draw directed edges (hybrid edges) with curvature
    directed_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for edge_tuple in network._graph.directed_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            edge_key = (u, v)
            if edge_key not in directed_groups:
                directed_groups[edge_key] = []
            directed_groups[edge_key].append((u, v, key))
    
    for (u, v), edges in directed_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=hybrid_edge_color, width=edge_width,
                arrows=True, arrowsize=arrow_size,
                arrowstyle='->', alpha=0.7, **kwargs
            )
        else:
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=hybrid_edge_color, linewidth=edge_width,
                    arrowstyle='->', mutation_scale=arrow_size,
                    alpha=0.7, zorder=2
                )
                ax.add_patch(arrow)
    
    # Draw nodes with different colors
    node_colors_map = {}
    node_sizes_map = {}
    for node in network._graph.nodes:
        if node in leaves:
            node_colors_map[node] = leaf_color
            node_sizes_map[node] = leaf_size
        elif node in hybrid_nodes:
            node_colors_map[node] = hybrid_color
            node_sizes_map[node] = node_size
        else:
            node_colors_map[node] = node_color
            node_sizes_map[node] = node_size
    
    nx.draw_networkx_nodes(
        G_nx, pos, ax=ax, node_color=[node_colors_map[n] for n in G_nx.nodes],
        node_size=[node_sizes_map[n] for n in G_nx.nodes], alpha=0.9, **kwargs
    )
    
    # Draw labels positioned on outside
    if with_labels:
        labels = {}
        for node in network._graph.nodes:
            label = network.get_label(node)
            if label:
                labels[node] = label
            else:
                labels[node] = str(node)
        
        # Position labels on outside (radial from center)
        label_pos = _position_labels_outside(pos, leaves)
        
        nx.draw_networkx_labels(G_nx, label_pos, labels=labels, ax=ax, font_size=10, font_weight='bold')
    
    network_type = 'Semi-Directed' if type(network).__name__ == 'SemiDirectedPhyNetwork' else 'Mixed'
    ax.set_title(f'{network_type} Phylogenetic Network', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return ax


def _hierarchical_layout_top_bottom(
    network: 'DirectedPhyNetwork',
    G_nx: nx.DiGraph
) -> dict[T, tuple[float, float]]:
    """
    Compute hierarchical layout for directed network (top to bottom).
    
    Root at top, leaves at bottom, minimizing crossings.
    """
    # Get root and compute topological levels
    root = network.root_node
    leaves = network.leaves
    
    # Compute levels using BFS from root
    levels: dict[T, int] = {}
    queue = [(root, 0)]
    visited = {root}
    
    while queue:
        node, level = queue.pop(0)
        levels[node] = level
        
        for child in network.children(node):
            if child not in visited:
                visited.add(child)
                queue.append((child, level + 1))
    
    # Group nodes by level
    level_groups: dict[int, list[T]] = {}
    for node, level in levels.items():
        if level not in level_groups:
            level_groups[level] = []
        level_groups[level].append(node)
    
    # Compute positions
    pos = {}
    max_level = max(levels.values()) if levels else 0
    y_spacing = 1.0 / (max_level + 1) if max_level > 0 else 1.0
    
    for level in sorted(level_groups.keys()):
        nodes_at_level = sorted(level_groups[level], key=lambda n: (network.get_label(n) or str(n)))
        num_nodes = len(nodes_at_level)
        x_spacing = 1.0 / (num_nodes + 1) if num_nodes > 0 else 0.5
        
        for idx, node in enumerate(nodes_at_level):
            x = (idx + 1) * x_spacing
            y = 1.0 - level * y_spacing
            pos[node] = (x, y)
    
    return pos


def _hierarchical_layout_left_right(
    network: 'DirectedPhyNetwork',
    G_nx: nx.DiGraph
) -> dict[T, tuple[float, float]]:
    """
    Compute hierarchical layout for directed network (left to right).
    
    Root at left, leaves at right, minimizing crossings.
    """
    # Get top-bottom layout and rotate 90 degrees
    pos_tb = _hierarchical_layout_top_bottom(network, G_nx)
    pos = {}
    for node, (x, y) in pos_tb.items():
        pos[node] = (y, 1.0 - x)  # Rotate 90 degrees clockwise
    
    return pos


def _position_labels_outside(
    pos: dict[T, tuple[float, float]],
    leaves: set[T]
) -> dict[T, tuple[float, float]]:
    """
    Position labels on outside of layout (radial from center).
    
    Parameters
    ----------
    pos : dict
        Node positions.
    leaves : set
        Set of leaf nodes (labels positioned further out).
    
    Returns
    -------
    dict
        Label positions.
    """
    if not pos:
        return {}
    
    # Compute center
    all_x = [x for x, y in pos.values()]
    all_y = [y for x, y in pos.values()]
    center_x = (min(all_x) + max(all_x)) / 2
    center_y = (min(all_y) + max(all_y)) / 2
    
    label_pos = {}
    offset_multiplier = 1.3  # Base offset
    leaf_offset_multiplier = 1.5  # Extra offset for leaves
    
    for node, (x, y) in pos.items():
        # Vector from center to node
        dx = x - center_x
        dy = y - center_y
        dist = np.sqrt(dx**2 + dy**2) if (dx != 0 or dy != 0) else 1.0
        
        if dist > 0:
            # Normalize
            dx_norm = dx / dist
            dy_norm = dy / dist
        else:
            dx_norm = 1.0
            dy_norm = 0.0
        
        # Offset distance
        if node in leaves:
            offset = dist * leaf_offset_multiplier * 0.15
        else:
            offset = dist * offset_multiplier * 0.1
        
        label_pos[node] = (x + dx_norm * offset, y + dy_norm * offset)
    
    return label_pos


def plot_tree(
    tree: 'DirectedPhyNetwork',
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """
    Plot a phylogenetic tree (DirectedPhyNetwork that is a tree).
    
    This is a convenience function that calls plot_network with tree-specific defaults.
    
    Parameters
    ----------
    tree : DirectedPhyNetwork
        The phylogenetic tree to plot (must be a tree, not a network).
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    **kwargs
        Additional keyword arguments passed to plot_network.
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    """
    return plot_network(tree, ax=ax, **kwargs)


# ========== New Layout-Based Plotting Functions ==========


def plot_network_with_layout(
    network: 'DirectedPhyNetwork | MixedPhyNetwork | SemiDirectedPhyNetwork',
    ax: plt.Axes | None = None,
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 500,
    leaf_size: int = 600,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    label_offset: float = 0.1,
    orientation: str = 'top-bottom',
    use_branch_lengths: bool = False,
    backend: str = 'matplotlib',
    **kwargs
) -> plt.Axes | Any:
    """
    Plot a phylogenetic network using optimal layout algorithms.
    
    This function uses the new layout computation functions for better
    crossing minimization and parallel edge handling. For DirectedPhyNetwork,
    it uses hierarchical layout with crossing minimization. For SemiDirectedPhyNetwork,
    it uses force-directed layout with leaf repulsion.
    
    Parameters
    ----------
    network : DirectedPhyNetwork | MixedPhyNetwork | SemiDirectedPhyNetwork
        The phylogenetic network to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, computes optimal layout.
    node_color : str, optional
        Color for internal nodes. Default is 'lightblue'.
    leaf_color : str, optional
        Color for leaf nodes. Default is 'lightgreen'.
    hybrid_color : str, optional
        Color for hybrid nodes. Default is 'salmon'.
    node_size : int, optional
        Size of internal nodes. Default is 500.
    leaf_size : int, optional
        Size of leaf nodes. Default is 600.
    edge_color : str, optional
        Color for regular edges. Default is 'gray'.
    hybrid_edge_color : str, optional
        Color for hybrid edges. Default is 'red'.
    edge_width : float, optional
        Width of edges. Default is 2.0.
    arrow_size : int, optional
        Size of arrows. Default is 20.
    with_labels : bool, optional
        Whether to show node labels. Default is True.
    label_offset : float, optional
        Offset for labels from nodes. Default is 0.1.
    orientation : str, optional
        Layout orientation: 'top-bottom' or 'left-right'. Default is 'top-bottom'.
    use_branch_lengths : bool, optional
        If True, scale node positions based on branch lengths (DirectedPhyNetwork only).
        Default is False.
    backend : str, optional
        Backend to use: 'matplotlib' or 'pyqtgraph'. Default is 'matplotlib'.
        If pyqtgraph is not available, falls back to matplotlib.
    **kwargs
        Additional keyword arguments passed to plotting functions.
    
    Returns
    -------
    matplotlib.axes.Axes | Any
        The axes object (matplotlib) or plot widget (pyqtgraph).
    
    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.network_plot import plot_network_with_layout
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> plot_network_with_layout(net)
    <Axes: >
    """
    # Auto-detect backend if pyqtgraph is available
    if backend == 'pyqtgraph' and not HAS_PYQTGAPH:
        backend = 'matplotlib'
    
    # Determine network type and call appropriate function
    network_type = type(network).__name__
    
    if network_type == 'DirectedPhyNetwork':
        if backend == 'pyqtgraph':
            return plot_directed_network_pyqtgraph(
                network, pos=pos, node_color=node_color, leaf_color=leaf_color,
                hybrid_color=hybrid_color, node_size=node_size, leaf_size=leaf_size,
                edge_color=edge_color, hybrid_edge_color=hybrid_edge_color,
                edge_width=edge_width, arrow_size=arrow_size, with_labels=with_labels,
                label_offset=label_offset, orientation=orientation,
                use_branch_lengths=use_branch_lengths, **kwargs
            )
        else:
            return plot_directed_network_with_layout(
                network, ax=ax, pos=pos, node_color=node_color, leaf_color=leaf_color,
                hybrid_color=hybrid_color, node_size=node_size, leaf_size=leaf_size,
                edge_color=edge_color, hybrid_edge_color=hybrid_edge_color,
                edge_width=edge_width, arrow_size=arrow_size, with_labels=with_labels,
                label_offset=label_offset, orientation=orientation,
                use_branch_lengths=use_branch_lengths, **kwargs
            )
    elif network_type in ('MixedPhyNetwork', 'SemiDirectedPhyNetwork'):
        if backend == 'pyqtgraph':
            return plot_semidirected_network_pyqtgraph(
                network, pos=pos, node_color=node_color, leaf_color=leaf_color,
                hybrid_color=hybrid_color, node_size=node_size, leaf_size=leaf_size,
                edge_color=edge_color, hybrid_edge_color=hybrid_edge_color,
                edge_width=edge_width, arrow_size=arrow_size, with_labels=with_labels,
                label_offset=label_offset, **kwargs
            )
        else:
            return plot_semidirected_network_with_layout(
                network, ax=ax, pos=pos, node_color=node_color, leaf_color=leaf_color,
                hybrid_color=hybrid_color, node_size=node_size, leaf_size=leaf_size,
                edge_color=edge_color, hybrid_edge_color=hybrid_edge_color,
                edge_width=edge_width, arrow_size=arrow_size, with_labels=with_labels,
                label_offset=label_offset, **kwargs
            )
    else:
        raise ValueError(f"Unsupported network type: {network_type}")


def plot_directed_network_with_layout(
    network: 'DirectedPhyNetwork',
    ax: plt.Axes | None = None,
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 500,
    leaf_size: int = 600,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    label_offset: float = 0.1,
    orientation: str = 'top-bottom',
    use_branch_lengths: bool = False,
    **kwargs
) -> plt.Axes:
    """
    Plot a DirectedPhyNetwork using optimal hierarchical layout.
    
    Uses `compute_directed_layout()` for crossing minimization and optimal positioning.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, computes optimal layout.
    orientation : str, optional
        'top-bottom' (default) or 'left-right'.
    use_branch_lengths : bool, optional
        If True, scale positions based on branch lengths. Default is False.
    **kwargs
        Additional keyword arguments (see plot_network_with_layout for full list).
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    
    # Get network properties
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Compute optimal layout if not provided
    if pos is None:
        pos = compute_directed_layout(
            network,
            orientation=orientation,
            use_branch_lengths=use_branch_lengths,
        )
    
    # Create NetworkX DiGraph for drawing
    G_nx = nx.DiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        G_nx.add_edge(u, v)
    
    # Group edges by (u, v) for parallel edge handling
    edge_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for u, v, key in network._graph.edges(keys=True):
        edge_key = (u, v)
        if edge_key not in edge_groups:
            edge_groups[edge_key] = []
        edge_groups[edge_key].append((u, v, key))
    
    # Draw edges with curvature for parallel edges
    for (u, v), edges in edge_groups.items():
        is_hybrid_edge = v in hybrid_nodes
        edge_col = hybrid_edge_color if is_hybrid_edge else edge_color
        
        num_parallel = len(edges)
        
        if num_parallel == 1:
            # Single edge - draw straight
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=edge_col, width=edge_width,
                arrows=True, arrowsize=arrow_size,
                arrowstyle='->', alpha=0.7, **kwargs
            )
        else:
            # Multiple parallel edges - draw with curvature
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                # Create curved path
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=edge_col, linewidth=edge_width,
                    arrowstyle='->', mutation_scale=arrow_size,
                    alpha=0.7, zorder=1
                )
                ax.add_patch(arrow)
    
    # Draw nodes with different colors for leaves and hybrids
    node_colors_map = {}
    node_sizes_map = {}
    for node in network._graph.nodes:
        if node in leaves:
            node_colors_map[node] = leaf_color
            node_sizes_map[node] = leaf_size
        elif node in hybrid_nodes:
            node_colors_map[node] = hybrid_color
            node_sizes_map[node] = node_size
        else:
            node_colors_map[node] = node_color
            node_sizes_map[node] = node_size
    
    nx.draw_networkx_nodes(
        G_nx, pos, ax=ax, node_color=[node_colors_map[n] for n in G_nx.nodes],
        node_size=[node_sizes_map[n] for n in G_nx.nodes], alpha=0.9, **kwargs
    )
    
    # Draw labels
    if with_labels:
        labels = {}
        for node in network._graph.nodes:
            label = network.get_label(node)
            if label:
                labels[node] = label
            else:
                labels[node] = str(node)
        
        # Position labels offset from nodes
        label_pos = {}
        if orientation == 'top-bottom':
            for node, (x, y) in pos.items():
                if node in leaves:
                    label_pos[node] = (x, y - label_offset)
                else:
                    label_pos[node] = (x, y + label_offset)
        else:  # left-right
            for node, (x, y) in pos.items():
                if node in leaves:
                    label_pos[node] = (x + label_offset, y)
                else:
                    label_pos[node] = (x - label_offset, y)
        
        nx.draw_networkx_labels(G_nx, label_pos, labels=labels, ax=ax, font_size=10, font_weight='bold')
    
    ax.set_title('Directed Phylogenetic Network', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return ax


def plot_semidirected_network_with_layout(
    network: 'MixedPhyNetwork | SemiDirectedPhyNetwork',
    ax: plt.Axes | None = None,
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 500,
    leaf_size: int = 600,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    label_offset: float = 0.15,
    **kwargs
) -> plt.Axes:
    """
    Plot a MixedPhyNetwork or SemiDirectedPhyNetwork using optimal force-directed layout.
    
    Uses `compute_semidirected_layout()` for leaf repulsion and optimal positioning.
    
    Parameters
    ----------
    network : MixedPhyNetwork | SemiDirectedPhyNetwork
        The mixed/semi-directed phylogenetic network to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, computes optimal layout.
    **kwargs
        Additional keyword arguments (see plot_network_with_layout for full list).
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))
    
    # Get network properties
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Compute optimal layout if not provided
    if pos is None:
        pos = compute_semidirected_layout(network)
    
    # Create NetworkX MultiGraph for drawing
    G_nx = nx.MultiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    # Add all edges for drawing
    for edge_tuple in network._graph.directed_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            G_nx.add_edge(u, v, key=key)
    for edge_tuple in network._graph.undirected_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            G_nx.add_edge(u, v, key=key)
    
    # Draw undirected edges with curvature for parallel edges
    undirected_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for edge_tuple in network._graph.undirected_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            edge_key = (min(u, v), max(u, v))
            if edge_key not in undirected_groups:
                undirected_groups[edge_key] = []
            undirected_groups[edge_key].append((u, v, key))
    
    for (u, v), edges in undirected_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=edge_color, width=edge_width,
                style='dashed', alpha=0.7, **kwargs
            )
        else:
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=edge_color, linewidth=edge_width,
                    linestyle='dashed', alpha=0.7, zorder=1,
                    arrowstyle='-'
                )
                ax.add_patch(arrow)
    
    # Draw directed edges (hybrid edges) with curvature
    directed_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for edge_tuple in network._graph.directed_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            edge_key = (u, v)
            if edge_key not in directed_groups:
                directed_groups[edge_key] = []
            directed_groups[edge_key].append((u, v, key))
    
    for (u, v), edges in directed_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=hybrid_edge_color, width=edge_width,
                arrows=True, arrowsize=arrow_size,
                arrowstyle='->', alpha=0.7, **kwargs
            )
        else:
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=hybrid_edge_color, linewidth=edge_width,
                    arrowstyle='->', mutation_scale=arrow_size,
                    alpha=0.7, zorder=2
                )
                ax.add_patch(arrow)
    
    # Draw nodes with different colors
    node_colors_map = {}
    node_sizes_map = {}
    for node in network._graph.nodes:
        if node in leaves:
            node_colors_map[node] = leaf_color
            node_sizes_map[node] = leaf_size
        elif node in hybrid_nodes:
            node_colors_map[node] = hybrid_color
            node_sizes_map[node] = node_size
        else:
            node_colors_map[node] = node_color
            node_sizes_map[node] = node_size
    
    nx.draw_networkx_nodes(
        G_nx, pos, ax=ax, node_color=[node_colors_map[n] for n in G_nx.nodes],
        node_size=[node_sizes_map[n] for n in G_nx.nodes], alpha=0.9, **kwargs
    )
    
    # Draw labels positioned on outside
    if with_labels:
        labels = {}
        for node in network._graph.nodes:
            label = network.get_label(node)
            if label:
                labels[node] = label
            else:
                labels[node] = str(node)
        
        # Position labels on outside (radial from center)
        label_pos = _position_labels_outside(pos, leaves)
        
        nx.draw_networkx_labels(G_nx, label_pos, labels=labels, ax=ax, font_size=10, font_weight='bold')
    
    network_type = 'Semi-Directed' if type(network).__name__ == 'SemiDirectedPhyNetwork' else 'Mixed'
    ax.set_title(f'{network_type} Phylogenetic Network', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return ax


# ========== PyQtGraph Backend (Optional) ==========


def plot_directed_network_pyqtgraph(
    network: 'DirectedPhyNetwork',
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 10,
    leaf_size: int = 12,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    with_labels: bool = True,
    orientation: str = 'top-bottom',
    use_branch_lengths: bool = False,
    **kwargs
) -> Any:
    """
    Plot a DirectedPhyNetwork using PyQtGraph backend.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to plot.
    pos : dict, optional
        Dictionary of node positions. If None, computes optimal layout.
    **kwargs
        Additional keyword arguments (see plot_network_with_layout for full list).
    
    Returns
    -------
    Any
        PyQtGraph plot widget.
    
    Raises
    ------
    ImportError
        If pyqtgraph is not installed.
    """
    if not HAS_PYQTGAPH:
        raise ImportError("pyqtgraph is not installed. Install it with: pip install pyqtgraph")
    
    # Compute optimal layout if not provided
    if pos is None:
        pos = compute_directed_layout(
            network,
            orientation=orientation,
            use_branch_lengths=use_branch_lengths,
        )
    
    # Create plot widget
    win = pg.GraphicsLayoutWidget(show=True, title="Directed Phylogenetic Network")
    plot = win.addPlot(title="Directed Phylogenetic Network")
    
    # Get network properties
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Extract positions
    x_coords = [pos[node][0] for node in network._graph.nodes]
    y_coords = [pos[node][1] for node in network._graph.nodes]
    
    # Plot nodes with different colors
    for node in network._graph.nodes:
        x, y = pos[node]
        if node in leaves:
            color = leaf_color
            size = leaf_size
        elif node in hybrid_nodes:
            color = hybrid_color
            size = node_size
        else:
            color = node_color
            size = node_size
        
        plot.plot([x], [y], pen=None, symbol='o', symbolSize=size, symbolBrush=color)
    
    # Plot edges
    for u, v, key in network._graph.edges(keys=True):
        is_hybrid_edge = v in hybrid_nodes
        edge_col = hybrid_edge_color if is_hybrid_edge else edge_color
        
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        plot.plot([x1, x2], [y1, y2], pen={'color': edge_col, 'width': edge_width})
    
    # Add labels if requested
    if with_labels:
        for node in network._graph.nodes:
            label = network.get_label(node) or str(node)
            x, y = pos[node]
            text = pg.TextItem(label)
            text.setPos(x, y)
            plot.addItem(text)
    
    return win


def plot_semidirected_network_pyqtgraph(
    network: 'MixedPhyNetwork | SemiDirectedPhyNetwork',
    pos: dict[T, tuple[float, float]] | None = None,
    node_color: str = 'lightblue',
    leaf_color: str = 'lightgreen',
    hybrid_color: str = 'salmon',
    node_size: int = 10,
    leaf_size: int = 12,
    edge_color: str = 'gray',
    hybrid_edge_color: str = 'red',
    edge_width: float = 2.0,
    with_labels: bool = True,
    **kwargs
) -> Any:
    """
    Plot a MixedPhyNetwork or SemiDirectedPhyNetwork using PyQtGraph backend.
    
    Parameters
    ----------
    network : MixedPhyNetwork | SemiDirectedPhyNetwork
        The mixed/semi-directed phylogenetic network to plot.
    pos : dict, optional
        Dictionary of node positions. If None, computes optimal layout.
    **kwargs
        Additional keyword arguments (see plot_network_with_layout for full list).
    
    Returns
    -------
    Any
        PyQtGraph plot widget.
    
    Raises
    ------
    ImportError
        If pyqtgraph is not installed.
    """
    if not HAS_PYQTGAPH:
        raise ImportError("pyqtgraph is not installed. Install it with: pip install pyqtgraph")
    
    # Compute optimal layout if not provided
    if pos is None:
        pos = compute_semidirected_layout(network)
    
    # Create plot widget
    network_type = 'Semi-Directed' if type(network).__name__ == 'SemiDirectedPhyNetwork' else 'Mixed'
    win = pg.GraphicsLayoutWidget(show=True, title=f"{network_type} Phylogenetic Network")
    plot = win.addPlot(title=f"{network_type} Phylogenetic Network")
    
    # Get network properties
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Plot nodes with different colors
    for node in network._graph.nodes:
        x, y = pos[node]
        if node in leaves:
            color = leaf_color
            size = leaf_size
        elif node in hybrid_nodes:
            color = hybrid_color
            size = node_size
        else:
            color = node_color
            size = node_size
        
        plot.plot([x], [y], pen=None, symbol='o', symbolSize=size, symbolBrush=color)
    
    # Plot undirected edges
    for edge_tuple in network._graph.undirected_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            plot.plot([x1, x2], [y1, y2], pen={'color': edge_color, 'width': edge_width, 'style': 2})  # Dashed
    
    # Plot directed edges (hybrid edges)
    for edge_tuple in network._graph.directed_edges_iter(keys=True):
        if len(edge_tuple) >= 3:
            u, v, key = edge_tuple[:3]
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            plot.plot([x1, x2], [y1, y2], pen={'color': hybrid_edge_color, 'width': edge_width})
    
    # Add labels if requested
    if with_labels:
        for node in network._graph.nodes:
            label = network.get_label(node) or str(node)
            x, y = pos[node]
            text = pg.TextItem(label)
            text.setPos(x, y)
            plot.addItem(text)
    
    return win
