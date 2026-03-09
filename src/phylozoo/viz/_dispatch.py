"""
Type dispatch for viz.plot().

Maps supported object types to their plotter functions and validates
layout support per type.
"""

from __future__ import annotations

from typing import Any, Callable

from phylozoo.utils.exceptions import PhyloZooLayoutError, PhyloZooTypeError

# Lazy imports to avoid circular deps
_DNET_PLOTTER: Callable[..., Any] | None = None
_SDNET_PLOTTER: Callable[..., Any] | None = None
_DMGRAPH_PLOTTER: Callable[..., Any] | None = None
_MMGRAPH_PLOTTER: Callable[..., Any] | None = None


def _get_plotter(obj: Any) -> tuple[Callable[..., Any], type[Any]]:
    """Return (plotter_fn, obj_type) for the given object."""
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

    global _DNET_PLOTTER, _SDNET_PLOTTER, _DMGRAPH_PLOTTER, _MMGRAPH_PLOTTER

    if _DNET_PLOTTER is None:
        from phylozoo.viz.dnetwork.plot import plot_dnetwork
        _DNET_PLOTTER = plot_dnetwork
    if _SDNET_PLOTTER is None:
        from phylozoo.viz.sdnetwork.plot import plot_sdnetwork
        _SDNET_PLOTTER = plot_sdnetwork
    if _DMGRAPH_PLOTTER is None:
        from phylozoo.viz.d_multigraph.plot import plot_dmgraph
        _DMGRAPH_PLOTTER = plot_dmgraph
    if _MMGRAPH_PLOTTER is None:
        from phylozoo.viz.m_multigraph.plot import plot_mmgraph
        _MMGRAPH_PLOTTER = plot_mmgraph

    if isinstance(obj, DirectedPhyNetwork):
        return _DNET_PLOTTER, DirectedPhyNetwork
    if isinstance(obj, SemiDirectedPhyNetwork):
        return _SDNET_PLOTTER, SemiDirectedPhyNetwork
    if isinstance(obj, DirectedMultiGraph):
        return _DMGRAPH_PLOTTER, DirectedMultiGraph
    if isinstance(obj, MixedMultiGraph):
        return _MMGRAPH_PLOTTER, MixedMultiGraph

    raise PhyloZooTypeError(
        f"Cannot plot {type(obj).__name__}. "
        "Supported types: DirectedPhyNetwork, SemiDirectedPhyNetwork, "
        "DirectedMultiGraph, MixedMultiGraph."
    )


def resolve_layout(obj: Any, layout: str) -> str:
    """
    Resolve layout='auto' to the default for the object type.
    Validate that the requested layout is supported for this type.
    """
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

    if layout != 'auto':
        resolved = layout
    elif isinstance(obj, DirectedPhyNetwork):
        resolved = 'pz-dag'
    elif isinstance(obj, SemiDirectedPhyNetwork):
        resolved = 'twopi'
    elif isinstance(obj, (DirectedMultiGraph, MixedMultiGraph)):
        resolved = 'spring'
    else:
        resolved = 'spring'

    # Validate layout support per type
    if isinstance(obj, DirectedPhyNetwork):
        if resolved.startswith('pz-') and resolved != 'pz-dag':
            raise PhyloZooLayoutError(
                f"Layout '{resolved}' is not supported for DirectedPhyNetwork. "
                "Use 'pz-dag' or a generic layout (spring, circular, dot, etc.)."
            )
    elif isinstance(obj, SemiDirectedPhyNetwork):
        if resolved.startswith('pz-') and resolved != 'pz-radial':
            raise PhyloZooLayoutError(
                f"Layout '{resolved}' is not supported for SemiDirectedPhyNetwork. "
                "Use 'pz-radial' (trees only) or a generic layout (twopi, spring, etc.)."
            )
    elif isinstance(obj, (DirectedMultiGraph, MixedMultiGraph)):
        if resolved.startswith('pz-'):
            raise PhyloZooLayoutError(
                f"PhyloZoo layouts (pz-*) are not supported for {type(obj).__name__}. "
                "Use a generic layout: spring, circular, dot, neato, twopi, etc."
            )

    return resolved
