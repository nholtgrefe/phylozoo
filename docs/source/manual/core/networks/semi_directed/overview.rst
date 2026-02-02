Overview
========

The :mod:`phylozoo.core.network.sdnetwork` module provides classes and functions for working with
semi-directed phylogenetic networks. Semi-directed networks are networks with directed hybrid edges
and undirected tree edges, allowing for flexible representation without requiring a fixed root.
Hybrid nodes have in-degree >= 2 and total degree = in-degree + 1. Semi-directed networks cannot have undirected cycles.

All classes and functions on this page can be imported from the core network module:

.. code-block:: python

   from phylozoo.core.network.sdnetwork import *
   # or directly
   from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

Documentation
-------------

* :doc:`Semi-Directed Network (Class) <semi_directed_network_class>`: The SemiDirectedPhyNetwork class and its properties
* :doc:`Semi-Directed Network (Advanced Features) <semi_directed_network_algorithms>`: Network features, transformations, classifications, isomorphism
* :doc:`Semi-Directed Generators <generators>`: Level-k network generators

.. seealso::
   For directed networks, see :doc:`Directed Networks <../directed/overview>`.
