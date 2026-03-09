"""
Public API for SemiDirectedPhyNetwork plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt

from phylozoo.utils.exceptions import PhyloZooLayoutError
from phylozoo.viz._layout_utils import compute_layout_center

from .style import SDNetStyle, default_style
from .layout import compute_nx_layout, compute_pz_radial_layout
from phylozoo.viz._render import render_layout

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork


def plot_sdnetwork(
    network: 'SemiDirectedPhyNetwork',
    layout: str = 'twopi',
    style: SDNetStyle | None = None,
    ax: Any | None = None,
    show: bool = False,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a SemiDirectedPhyNetwork.

    This is the main public API function for plotting semi-directed networks.
    It handles layout computation, styling, and rendering using matplotlib.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The network to plot.
    layout : str, optional
        Layout algorithm. PhyloZoo: 'pz-radial' (trees only). NetworkX: 'spring',
        'circular', 'kamada_kawai', 'planar', 'random', 'shell', 'spectral', 'spiral',
        'bipartite'. Graphviz: 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'twopi'.
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
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz import plot
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
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

    if layout == 'pz-radial':
        computed_layout = compute_pz_radial_layout(network, **layout_kwargs)
    elif layout.startswith('pz-'):
        raise PhyloZooLayoutError(
            f"Unknown PhyloZoo layout: '{layout}'. "
            "Supported PhyloZoo layouts: 'pz-radial'"
        )
    else:
        computed_layout = compute_nx_layout(network, layout=layout, **layout_kwargs)

    network_obj = computed_layout.network
    positions = computed_layout.positions
    center = compute_layout_center(positions)
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes

    def get_node_type(node: Any) -> str:
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
        radial_labels_for_leaves=(layout == 'pz-radial'),
    )

    if show:
        plt.show()

    return ax
