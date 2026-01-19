Splits
======

The split module provides classes for working with phylogenetic splits and split systems.
A split is a bipartition of a set of taxa, representing a division of the taxa into two
non-empty subsets. Split systems are collections of splits that can represent phylogenetic
trees or networks.

Main Classes
------------

A split is a 2-partition :math:`\{A, B\}` of a set of elements where :math:`A \cup B` equals
the full set and :math:`A \cap B = \emptyset`.

.. automodule:: phylozoo.core.split.base
   :members:
   :show-inheritance:

Split Systems
-------------

A split system is a collection of splits where each split covers the complete set of elements.
Weighted split systems assign positive weights to each split.

.. automodule:: phylozoo.core.split.splitsystem
   :members:
   :show-inheritance:

.. automodule:: phylozoo.core.split.weighted_splitsystem
   :members:
   :show-inheritance:

Algorithms
----------

Algorithms for working with split systems, including conversion to phylogenetic networks
and computation of distance matrices.

.. automodule:: phylozoo.core.split.algorithms
   :members:
   :show-inheritance:

Classification Functions
------------------------

Functions for classifying split systems, such as checking pairwise compatibility.

.. automodule:: phylozoo.core.split.classifications
   :members:
   :show-inheritance:

I/O Support
-----------

Split systems support reading and writing in NEXUS format.

.. automodule:: phylozoo.core.split.io
   :members:
   :show-inheritance:
