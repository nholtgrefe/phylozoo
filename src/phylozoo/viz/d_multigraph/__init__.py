"""DirectedMultiGraph plotting module."""

from .layout import DMGraphLayout, compute_nx_layout
from .style import DMGraphStyle, default_style

__all__ = [
    "compute_nx_layout",
    "plot_dmgraph",
    "DMGraphLayout",
    "DMGraphStyle",
    "default_style",
]


def __getattr__(name: str) -> object:
    """Import the matplotlib-based plotter only when it is asked for."""
    if name == "plot_dmgraph":
        from .plot import plot_dmgraph

        return plot_dmgraph
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
