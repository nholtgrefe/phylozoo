"""
Public API for DirectedPhyNetwork plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from phylozoo.utils.exceptions import PhyloZooLayoutError
from phylozoo.viz._matplotlib import plt

from phylozoo.viz._layout_utils import compute_layout_center
from .style import DNetStyle, default_style
from .layout import (
    compute_nx_layout,
    compute_pz_unrooted_layout,
    compute_pz_cladogram_layout,
    compute_pz_layered_layout,
    compute_pz_radial_layout,
)
from phylozoo.viz._render import render_layout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork


def plot_dnetwork(
    network: "DirectedPhyNetwork",
    layout: str = "pz-cladogram",
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
        Layout algorithm. PhyloZoo: 'pz-cladogram' (layered with a tree backbone;
        ``rectangular=False`` for straight edges), 'pz-layered' (layered, dot-like),
        'pz-radial' (circular cladogram) or 'pz-unrooted' (unrooted tree of blobs).
        NetworkX: 'spring', 'circular',
        'kamada_kawai', 'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite'.
        Graphviz: 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'pz-cladogram'.
    style : DNetStyle, optional
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
    if layout == "pz-unrooted" and (style.root_color is None or style.root_size is None):
        # In the unrooted drawing nothing else marks the root: make it stand out.
        style = style.copy()
        style.root_color = style.root_color or "#f5c542"
        style.root_size = style.root_size or 2.2 * style.node_size

    if ax is None:
        _, ax = plt.subplots()

    if layout == "pz-cladogram":
        computed_layout = compute_pz_cladogram_layout(network, **layout_kwargs)
    elif layout == "pz-layered":
        computed_layout = compute_pz_layered_layout(network, **layout_kwargs)
    elif layout == "pz-radial":
        computed_layout = compute_pz_radial_layout(network, **layout_kwargs)
    elif layout == "pz-unrooted":
        computed_layout = compute_pz_unrooted_layout(network, **layout_kwargs)
    elif layout.startswith("pz-"):
        raise PhyloZooLayoutError(
            f"Unknown PhyloZoo layout: '{layout}'. Supported PhyloZoo layouts: "
            "'pz-cladogram', 'pz-layered', 'pz-radial', 'pz-unrooted'"
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
            return "root"
        if node in leaves:
            return "leaf"
        if node in hybrid_nodes:
            return "hybrid"
        return "tree"

    layered = layout in ("pz-cladogram", "pz-layered")
    arrows = style.arrows or ("hybrid" if layered or layout == "pz-radial" else "all")
    leaf_label_direction: tuple[float, float] | None = None
    if layered and layout_kwargs.get("align_leaves", layout != "pz-layered"):
        # Leaves form the bottom (or right) layer: point their labels away from it.
        left_right = str(layout_kwargs.get("direction", "TD")).upper() == "LR"
        leaf_label_direction = (1.0, 0.0) if left_right else (0.0, -1.0)

    render_layout(
        ax,
        computed_layout.edge_routes,
        positions,
        style,
        center,
        get_node_type,
        network_obj.get_label,
        radial_labels_for_leaves=(layout == "pz-radial"),
        arrows=arrows,
        leaf_label_direction=leaf_label_direction,
    )

    if show:
        plt.show()

    return ax
