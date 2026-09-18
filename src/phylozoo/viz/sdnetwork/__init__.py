"""Plotting and styling for semi-directed phylogenetic networks (SemiDirectedPhyNetwork)."""

from .layout import (
    SDNetLayout,
    compute_nx_layout,
    compute_pz_unrooted_layout,
    compute_pz_radial_layout,
)
from .style import SDNetStyle, default_style

__all__ = [
    "compute_nx_layout",
    "compute_pz_unrooted_layout",
    "compute_pz_radial_layout",
    "plot_sdnetwork",
    "SDNetLayout",
    "SDNetStyle",
    "default_style",
]


def __getattr__(name: str) -> object:
    """Import the matplotlib-based plotter only when it is asked for."""
    if name == "plot_sdnetwork":
        from .plot import plot_sdnetwork

        return plot_sdnetwork
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
