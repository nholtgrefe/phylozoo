Overview
========

The :mod:`phylozoo.core.network` module provides classes for working with phylogenetic networks,
which represent evolutionary relationships including reticulate events such as hybridization,
recombination, and horizontal gene transfer.

PhyloZoo provides two main network types: **DirectedPhyNetwork** for rooted networks with all
edges directed, and **SemiDirectedPhyNetwork** for networks with directed hybrid edges and
undirected tree edges, allowing flexible representation without requiring a fixed root.


Submodules
----------

- :doc:`Directed Networks <directed/overview>` - Rooted networks with all edges directed
- :doc:`Semi-Directed Networks <semi_directed/overview>` - Networks with mixed directed and undirected edges

See Also
--------

- :doc:`API Reference <../../../../api/core/network/index>` - Complete function signatures and detailed examples
- :doc:`I/O <../../utils/io/index>` - File I/O operations and formats
- :doc:`Visualization <../../visualization/overview>` - Network visualization and plotting
