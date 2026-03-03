"""
Network visualization and plotting (viz).

The visualization module provides a flexible plotting system for phylogenetic
networks and graphs with support for multiple layout algorithms and customizable
styling. It uses a modular design with clear separation between layout
computation and rendering.

**When to use which plot function:**

- :func:`plot_dnetwork` – For :class:`DirectedPhyNetwork` (rooted phylogenetic
  networks). Use this for standard phylogenetic network visualization with
  taxon labels, hybrid nodes, and tree structure. Default layout: ``'pz-dag'``.

- :func:`plot_sdnetwork` – For :class:`SemiDirectedPhyNetwork` (semi-directed
  phylogenetic networks). Use for networks with undirected tree edges and
  directed hybrid edges. Default layout: ``'twopi'``.

- :func:`plot_dmgraph` – For :class:`DirectedMultiGraph` (low-level directed
  multigraph). Use when working with raw graph structures, e.g. when implementing
  algorithms or debugging. No semantic node types (root, leaf, hybrid).

- :func:`plot_mmgraph` – For :class:`MixedMultiGraph` (low-level mixed graph).
  Same as above but for graphs with both directed and undirected edges.
"""

from .dnetwork import plot_dnetwork
from .graphs import plot_dmgraph, plot_mmgraph
from .sdnetwork import plot_sdnetwork

__all__ = [
    'plot_dnetwork',
    'plot_sdnetwork',
    'plot_dmgraph',
    'plot_mmgraph',
]
