"""MixedMultiGraph plotting module."""

from .layout import MGraphLayout, compute_nx_layout
from .style import MGraphStyle, default_style

__all__ = [
    "compute_nx_layout",
    "plot_mmgraph",
    "MGraphLayout",
    "MGraphStyle",
    "default_style",
]


def __getattr__(name: str) -> object:
    """Import the matplotlib-based plotter only when it is asked for."""
    if name == "plot_mmgraph":
        from .plot import plot_mmgraph

        return plot_mmgraph
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
