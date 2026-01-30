Primitives
==========

Fundamental data structures used throughout PhyloZoo. These classes provide the underlying
graph and set-theoretic structures that support the higher-level phylogenetic network classes.

Partition
---------

A partition is a division of a set into non-empty, disjoint subsets (parts). Partitions are
used to represent splits and other set-theoretic structures in phylogenetic analysis.

.. automodule:: phylozoo.core.primitives.partition
   :members:
   :show-inheritance:

Circular Orderings
------------------

A circular ordering represents a cyclic arrangement of elements, such as the order of taxa
around a circular phylogenetic tree or network. Circular orderings are used in network
reconstruction algorithms and visualization.

.. automodule:: phylozoo.core.primitives.circular_ordering
   :members:
   :show-inheritance:
