Overview
========

The :mod:`phylozoo.core` module provides fundamental data structures and classes used 
throughout PhyloZoo for representing and analyzing phylogenetic networks, quartets, splits,
sequences, distance matrices, and other core phylogenetic concepts.

All core data structures are available at the package top-level for convenience, but functions must be imported from the submodules.

.. code-block:: python

   from phylozoo import (
       DirectedPhyNetwork,
       SemiDirectedPhyNetwork,
       Quartet,
       QuartetProfile,
       QuartetProfileSet,
       Split,
       SplitSystem,
       DistanceMatrix,
       MSA,
   )

Submodules
----------

The :mod:`phylozoo.core` module is organized into several submodules, each focusing on a 
specific data type or functionality:

- :doc:`Networks <networks/overview>` - Directed and semi-directed phylogenetic network classes
- :doc:`Quartets <quartets/overview>` - Quartet information
- :doc:`Splits <splits/overview>` - Bipartitions of taxa and split systems
- :doc:`Sequences <sequences>` - Multiple sequence alignment (MSA) handling
- :doc:`Distance Matrices <distance>` - Distance matrix classes and operations
- :doc:`Primitives <primitives/overview>` - Fundamental data structures (Partition, CircularOrdering, MixedMultiGraph,etc.)

These components are designed to be interoperable, allowing for seamless transitions 
between different data representations.

See Also
--------

- :doc:`API Reference <../../api/core/index>` — Complete function signatures and detailed examples
