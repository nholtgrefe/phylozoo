"""
Unified plotting API for phylogenetic networks and graphs.

This module provides a single plot() function that dispatches by object type.
"""

from __future__ import annotations

from typing import Any

from ._dispatch import _get_plotter, resolve_layout


def plot(
    obj: Any,
    layout: str = 'auto',
    style: Any = None,
    ax: Any | None = None,
    show: bool = False,
    **kwargs: Any,
) -> Any:
    """
    Plot a phylogenetic network or graph.

    Dispatches by object type: DirectedPhyNetwork, SemiDirectedPhyNetwork,
    DirectedMultiGraph, or MixedMultiGraph. Use layout='auto' for the default
    layout per type, or specify a layout name.

    Parameters
    ----------
    obj : DirectedPhyNetwork | SemiDirectedPhyNetwork | DirectedMultiGraph | MixedMultiGraph
        The network or graph to plot.
    layout : str, optional
        Layout algorithm. Use 'auto' for the default per type.
        DirectedPhyNetwork: 'pz-dag' (default), or NetworkX ('spring', 'circular',
        'kamada_kawai', 'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite')
        or Graphviz ('dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo').
        SemiDirectedPhyNetwork: 'twopi' (default), 'pz-radial' (trees only), or NetworkX/Graphviz.
        DirectedMultiGraph / MixedMultiGraph: NetworkX or Graphviz only (no PhyloZoo layouts).
        By default 'auto'.
    style : DNetStyle | SDNetStyle | DMGraphStyle | MGraphStyle, optional
        Styling configuration. If None, uses default style for the object type.
        By default None.
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on. If None, creates new figure and axes.
        By default None.
    show : bool, optional
        If True, automatically display the plot using plt.show().
        By default False.
    **kwargs
        Additional parameters for layout computation (e.g. layer_gap, leaf_gap for pz-dag).

    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.

    Raises
    ------
    PhyloZooTypeError
        If obj is not a supported type.
    PhyloZooLayoutError
        If layout is not supported for the object type.

    Examples
    --------
    >>> from phylozoo import DirectedPhyNetwork
    >>> from phylozoo.viz import plot
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> ax = plot(net)
    >>> ax = plot(net, layout='spring', show=True)
    """
    plotter, _ = _get_plotter(obj)
    resolved_layout = resolve_layout(obj, layout)
    return plotter(obj, layout=resolved_layout, style=style, ax=ax, show=show, **kwargs)
