Overview
========

The :mod:`phylozoo.core.primitives` module provides fundamental data structures
that serve as building blocks for phylogenetic analysis in PhyloZoo. These classes
represent core mathematical concepts used throughout the library, including partitions,
circular orderings, and graph primitives that form the foundation for network-based
phylogenetics.

Some classes on this page can be imported from the core primitives module:

.. code-block:: python

   from phylozoo.core.primitives import *
   # or directly
   from phylozoo.core.primitives import (
       Partition,
       CircularOrdering,
       CircularSetOrdering,
   )

Graph classes and functions are imported from their submodules:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

Submodules
----------

- :doc:`Partition <partition>` - Set partitions and refinement relationships
- :doc:`Circular Ordering <circular_ordering>` - Circular arrangements of elements and sets
- :doc:`Directed Multi-Graph <directed_multigraph>` - Directed graphs with parallel edges
- :doc:`Mixed Multi-Graph <mixed_multigraph>` - Mixed graphs with both directed and undirected edges

See Also
--------

- :doc:`API Reference <../../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Networks <../networks/overview>` - Network classes that use these primitives
- :doc:`Splits <../splits/overview>` - Split systems based on partitions
