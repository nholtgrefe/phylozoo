"""
Public API for DirectedPhyNetwork plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt

from phylozoo.utils.exceptions import PhyloZooLayoutError
from phylozoo.viz._layout_utils import compute_layout_center

from .style import DNetStyle, default_style
from .layout import compute_nx_layout, compute_pz_dag_layout
from phylozoo.viz._render import render_layout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork


def plot_dnetwork(
    network: 'DirectedPhyNetwork',
    layout: str = 'pz-dag',
    style: DNetStyle | None = None,
    ax: Any | None = None,
    show: bool = False,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a DirectedPhyNetwork.

    This is the main public API function for plotting networks. It handles
    layout computation, styling, and rendering using matplotlib.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to plot.
    layout : str, optional
        Layout algorithm. PhyloZoo: 'pz-dag'. NetworkX: 'spring', 'circular',
        'kamada_kawai', 'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite'.
        Graphviz: 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'pz-dag'.
    style : NetworkStyle, optional
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
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz import plot
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> ax = plot(net)
    """
    if style is None:
        style = default_style()

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if layout == 'pz-dag':
        computed_layout = compute_pz_dag_layout(network, **layout_kwargs)
    elif layout.startswith('pz-'):
        raise PhyloZooLayoutError(
            f"Unknown PhyloZoo layout: '{layout}'. "
            "Supported PhyloZoo layouts: 'pz-dag'"
        )
    else:
        computed_layout = compute_nx_layout(network, layout=layout, **layout_kwargs)

    network_obj = computed_layout.network
    positions = computed_layout.positions
    center = compute_layout_center(positions)
    leaves = network_obj.leaves
    hybrid_nodes = network_obj.hybrid_nodes
    root = network_obj.root_node

    def get_node_type(node: Any) -> str:
        if node == root:
            return 'root'
        if node in leaves:
            return 'leaf'
        if node in hybrid_nodes:
            return 'hybrid'
        return 'tree'

    render_layout(
        ax,
        computed_layout.edge_routes,
        positions,
        style,
        center,
        get_node_type,
        network_obj.get_label,
        radial_labels_for_leaves=False,
    )

    if show:
        plt.show()

    return ax
