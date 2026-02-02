Overview
========

The :mod:`phylozoo.core.network.dnetwork` module provides classes and functions for working with
directed phylogenetic networks. Directed networks are rooted networks where all edges are directed,
explicitly representing the direction of evolutionary time. Note that these are often also referred
to as rooted phylogenetic networks. A rooted phylogenetic tree is a special case of a directed phylogenetic network.

All classes and functions on this page can be imported from the dnetwork module:

.. code-block:: python

   from phylozoo.core.network.dnetwork import *
   # or directly
   from phylozoo.core.network.dnetwork import DirectedPhyNetwork

.. tip::
   For convenience, you can also use the alias ``phylozoo.core.dnetwork`` instead of
   ``phylozoo.core.network.dnetwork``. Both import paths work identically.

Documentation
-------------

- :doc:`Directed Network (Class) <directed_network_class>` - The DirectedPhyNetwork class and its properties
- :doc:`Directed Network (Advanced Features) <directed_network_algorithms>` - Network features, transformations, classifications, isomorphism
- :doc:`Directed Generator <directed_generator>` - Level-k network generators

See Also
--------

- :doc:`API Reference <../../../api/core/network>` - Complete function signatures and detailed examples
- :doc:`Semi-Directed Networks <../semi_directed/overview>` - Networks with mixed directed and undirected edges 
