Overview
========

The :mod:`phylozoo.core.network.dnetwork` module provides classes and functions for working with
directed phylogenetic networks. Directed networks are networks where all edges are directed,
representing the direction of evolutionary time. Hybrid nodes have in-degree >= 2 and out-degree 1.

All classes and functions on this page can be imported from the core network module:

.. code-block:: python

   from phylozoo.core.network.dnetwork import *
   # or directly
   from phylozoo.core.network.dnetwork import DirectedPhyNetwork

Documentation
-------------

* :doc:`Directed Network (Class) <directed_network_class>`: The DirectedPhyNetwork class and its properties
* :doc:`Directed Network (Advanced Features) <directed_network_algorithms>`: Network features, transformations, classifications, isomorphism
* :doc:`Directed Generators <directed_generator>`: Level-k network generators

.. seealso::
   For semi-directed networks, see :doc:`Semi-Directed Networks <../semi_directed/overview>`. 
