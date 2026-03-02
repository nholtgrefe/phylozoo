Overview
========

The :mod:`phylozoo.core.network.sdnetwork` module provides classes and functions for working with
semi-directed phylogenetic networks. Semi-directed networks have directed hybrid edges and undirected
tree edges, allowing for flexible representation without requiring a fixed root. Note that unrooted phylogenetic trees are a special case of a semi-directed phylogenetic network.

All classes and functions on this page can be imported from the sdnetwork module:

.. code-block:: python

   from phylozoo.core.network.sdnetwork import *
   # or directly
   from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

.. tip::
   For convenience, you can also use the alias ``phylozoo.core.sdnetwork`` instead of
   ``phylozoo.core.network.sdnetwork``. Both import paths work identically.

Documentation
-------------

- :doc:`Semi-Directed Network (Class) <semi_directed_network_class>` - The SemiDirectedPhyNetwork class and its properties
- :doc:`Semi-Directed Network (Advanced Features) <semi_directed_network_algorithms>` - Network features, transformations, classifications, isomorphism
- :doc:`Semi-Directed Generator <generators>` - Level-k network generators

See Also
--------

- :doc:`API Reference <../../../../api/core/network/index>` - Complete function signatures and detailed examples
- :doc:`Directed Networks <../directed/overview>` - Rooted networks with all edges directed
