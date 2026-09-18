"""Plotting and styling for directed phylogenetic networks (DirectedPhyNetwork)."""

from .layout import (
    DNetLayout,
    compute_nx_layout,
    compute_pz_unrooted_layout,
    compute_pz_cladogram_layout,
    compute_pz_layered_layout,
    compute_pz_radial_layout,
)
from .style import DNetStyle, default_style

__all__ = [
    "compute_nx_layout",
    "compute_pz_unrooted_layout",
    "compute_pz_cladogram_layout",
    "compute_pz_layered_layout",
    "compute_pz_radial_layout",
    "plot_dnetwork",
    "DNetLayout",
    "DNetStyle",
    "default_style",
]


def __getattr__(name: str) -> object:
    """Import the matplotlib-based plotter only when it is asked for."""
    if name == "plot_dnetwork":
        from .plot import plot_dnetwork

        return plot_dnetwork
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
