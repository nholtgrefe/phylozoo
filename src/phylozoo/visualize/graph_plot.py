"""
Graph plotting functions for DirectedMultiGraph and MixedMultiGraph.

This module provides functions for visualizing graph structures with support
for parallel edges displayed as curved arcs.
"""

from typing import Tuple, TypeVar
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import networkx as nx

T = TypeVar('T')


def plot_directed_multigraph(
    graph: 'DirectedMultiGraph',
    ax: plt.Axes | None = None,
    pos: dict | None = None,
    node_color: str = 'lightblue',
    node_size: int = 500,
    edge_color: str = 'blue',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    **kwargs
) -> plt.Axes:
    """
    Plot a DirectedMultiGraph with curved parallel edges.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multigraph to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, uses spring layout.
    node_color : str, optional
        Color for nodes. Default is 'lightblue'.
    node_size : int, optional
        Size of nodes. Default is 500.
    edge_color : str, optional
        Color for edges. Default is 'blue'.
    edge_width : float, optional
        Width of edges. Default is 2.0.
    arrow_size : int, optional
        Size of arrows. Default is 20.
    with_labels : bool, optional
        Whether to show node labels. Default is True.
    **kwargs
        Additional keyword arguments passed to networkx drawing functions.
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.visualize.graph_plot import plot_directed_multigraph
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2)
    0
    >>> G.add_edge(1, 2)  # Parallel edge
    1
    >>> plot_directed_multigraph(G)
    <Axes: >
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create NetworkX MultiDiGraph for layout
    G_nx = nx.MultiDiGraph()
    for node in graph.nodes():
        G_nx.add_node(node)
    for u, v, key in graph.edges(keys=True):
        G_nx.add_edge(u, v, key=key)
    
    # Compute layout if not provided
    if pos is None:
        pos = nx.spring_layout(G_nx, seed=42)
    
    # Draw nodes
    nx.draw_networkx_nodes(
        G_nx, pos, ax=ax, node_color=node_color,
        node_size=node_size, alpha=0.9, **kwargs
    )
    
    # Group edges by (u, v) to handle parallel edges
    edge_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for u, v, key in graph.edges(keys=True):
        edge_key = (u, v)
        if edge_key not in edge_groups:
            edge_groups[edge_key] = []
        edge_groups[edge_key].append((u, v, key))
    
    # Draw edges with curvature for parallel edges
    for (u, v), edges in edge_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            # Single edge - draw straight
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=edge_color, width=edge_width,
                arrows=True, arrowsize=arrow_size,
                arrowstyle='->', alpha=0.7, **kwargs
            )
        else:
            # Multiple parallel edges - draw with curvature
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                # Calculate curvature offset
                # Use different curvature for each parallel edge
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                
                # Get positions
                x1, y1 = pos[u]
                x2, y2 = pos[v]
                
                # Calculate control point for curved edge
                # Perpendicular to the line between nodes
                dx = x2 - x1
                dy = y2 - y1
                length = math.sqrt(dx**2 + dy**2)
                
                if length > 0:
                    # Perpendicular vector
                    perp_x = -dy / length
                    perp_y = dx / length
                    
                    # Control point offset
                    offset = curvature * length * 0.5
                    mid_x = (x1 + x2) / 2 + perp_x * offset
                    mid_y = (y1 + y2) / 2 + perp_y * offset
                    
                # Create curved path using FancyArrowPatch
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=edge_color, linewidth=edge_width,
                    arrowstyle='->', mutation_scale=arrow_size,
                    alpha=0.7, zorder=1
                )
                ax.add_patch(arrow)
    
    # Draw labels
    if with_labels:
        nx.draw_networkx_labels(G_nx, pos, ax=ax, font_size=10, font_weight='bold')
    
    ax.set_title('Directed MultiGraph', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return ax


def plot_mixed_multigraph(
    graph: 'MixedMultiGraph',
    ax: plt.Axes | None = None,
    pos: dict | None = None,
    node_color: str = 'lightblue',
    node_size: int = 500,
    directed_edge_color: str = 'red',
    undirected_edge_color: str = 'gray',
    edge_width: float = 2.0,
    arrow_size: int = 20,
    with_labels: bool = True,
    **kwargs
) -> plt.Axes:
    """
    Plot a MixedMultiGraph with curved parallel edges.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multigraph to plot.
    ax : matplotlib.axes.Axes, optional
        Matplotlib axes to plot on. If None, creates a new figure.
    pos : dict, optional
        Dictionary of node positions. If None, uses spring layout.
    node_color : str, optional
        Color for nodes. Default is 'lightblue'.
    node_size : int, optional
        Size of nodes. Default is 500.
    directed_edge_color : str, optional
        Color for directed edges. Default is 'red'.
    undirected_edge_color : str, optional
        Color for undirected edges. Default is 'gray'.
    edge_width : float, optional
        Width of edges. Default is 2.0.
    arrow_size : int, optional
        Size of arrows. Default is 20.
    with_labels : bool, optional
        Whether to show node labels. Default is True.
    **kwargs
        Additional keyword arguments passed to networkx drawing functions.
    
    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.visualize.graph_plot import plot_mixed_multigraph
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> G.add_undirected_edge(2, 3)  # Parallel undirected edge
    1
    >>> plot_mixed_multigraph(G)
    <Axes: >
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create NetworkX MultiGraph for layout (combining both directed and undirected)
    G_nx = nx.MultiGraph()
    for node in graph.nodes():
        G_nx.add_node(node)
    # Add all edges for layout calculation
    for u, v, key in graph._directed.edges(keys=True):
        G_nx.add_edge(u, v, key=key)
    for u, v, key in graph._undirected.edges(keys=True):
        G_nx.add_edge(u, v, key=key)
    
    # Compute layout if not provided
    if pos is None:
        pos = nx.spring_layout(G_nx, seed=42)
    
    # Draw nodes
    nx.draw_networkx_nodes(
        G_nx, pos, ax=ax, node_color=node_color,
        node_size=node_size, alpha=0.9, **kwargs
    )
    
    # Draw undirected edges with curvature for parallel edges
    undirected_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for u, v, key in graph._undirected.edges(keys=True):
        edge_key = (min(u, v), max(u, v))  # Canonical form
        if edge_key not in undirected_groups:
            undirected_groups[edge_key] = []
        undirected_groups[edge_key].append((u, v, key))
    
    for (u, v), edges in undirected_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            # Single edge - draw straight
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=undirected_edge_color, width=edge_width,
                style='dashed', alpha=0.7, **kwargs
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
                    color=undirected_edge_color, linewidth=edge_width,
                    linestyle='dashed', alpha=0.7, zorder=1,
                    arrowstyle='-'
                )
                ax.add_patch(arrow)
    
    # Draw directed edges with curvature for parallel edges
    directed_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for u, v, key in graph._directed.edges(keys=True):
        edge_key = (u, v)
        if edge_key not in directed_groups:
            directed_groups[edge_key] = []
        directed_groups[edge_key].append((u, v, key))
    
    for (u, v), edges in directed_groups.items():
        num_parallel = len(edges)
        
        if num_parallel == 1:
            # Single edge - draw straight
            nx.draw_networkx_edges(
                G_nx, pos, edgelist=[(u, v)], ax=ax,
                edge_color=directed_edge_color, width=edge_width,
                arrows=True, arrowsize=arrow_size,
                arrowstyle='->', alpha=0.7, **kwargs
            )
        else:
            # Multiple parallel edges - draw with curvature
            for idx, (u_edge, v_edge, key) in enumerate(edges):
                curvature = (idx - (num_parallel - 1) / 2) * 0.3
                
                x1, y1 = pos[u_edge]
                x2, y2 = pos[v_edge]
                
                # Create curved path with arrow
                arrow = mpatches.FancyArrowPatch(
                    (x1, y1), (x2, y2),
                    connectionstyle=f"arc3,rad={curvature}",
                    color=directed_edge_color, linewidth=edge_width,
                    arrowstyle='->', mutation_scale=arrow_size,
                    alpha=0.7, zorder=2
                )
                ax.add_patch(arrow)
    
    # Draw labels
    if with_labels:
        nx.draw_networkx_labels(G_nx, pos, ax=ax, font_size=10, font_weight='bold')
    
    ax.set_title('Mixed MultiGraph', fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return ax

