"""
Public API for MixedMultiGraph plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt

from phylozoo.utils.exceptions import PhyloZooLayoutError
from phylozoo.viz._layout_utils import compute_layout_center

from .layout.nx import compute_nx_layout
from .style import MGraphStyle, default_style
from phylozoo.viz._render import render_layout

if TYPE_CHECKING:
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph


def plot_mmgraph(
    graph: 'MixedMultiGraph',
    layout: str = 'spring',
    style: MGraphStyle | None = None,
    ax: Any | None = None,
    show: bool = False,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a MixedMultiGraph.

    This is the main public API function for plotting mixed multigraphs.
    It handles layout computation, styling, and rendering using matplotlib.

    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to plot.
    layout : str, optional
        Layout algorithm. NetworkX: 'spring', 'circular', 'kamada_kawai', 'planar',
        'random', 'shell', 'spectral', 'spiral', 'bipartite'. Graphviz: 'dot',
        'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'spring'.
    style : MGraphStyle, optional
        Styling configuration. If None, uses default style.
        By default None.
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on. If None, creates new figure and axes.
        By default None.
    show : bool, optional
        If True, automatically display the plot using plt.show().
        By default False.
    **layout_kwargs
        Additional parameters for layout computation.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported.

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.viz import plot
    >>>
    >>> G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])
    >>> ax = plot(G, layout='spring')
    """
    if style is None:
        style = default_style()

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    computed_layout = compute_nx_layout(graph, layout=layout, **layout_kwargs)
    positions = computed_layout.positions
    center = compute_layout_center(positions)

    def get_label(node: Any) -> str | None:
        if node in graph._directed.nodes:
            return graph._directed.nodes[node].get('label', str(node))
        if node in graph._undirected.nodes:
            return graph._undirected.nodes[node].get('label', str(node))
        return str(node)

    render_layout(
        ax,
        computed_layout.edge_routes,
        positions,
        style,
        center,
        get_node_type=lambda _: 'generic',
        get_label=get_label,
        radial_labels_for_leaves=False,
    )

    if show:
        plt.show()

    return ax
