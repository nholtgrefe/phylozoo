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

Directed Multi-Graph
--------------------

A directed multi-graph allows multiple directed edges between the same pair of vertices.
This is the underlying graph structure for directed phylogenetic networks, where multiple
edges can represent parallel arcs.

.. automodule:: phylozoo.core.primitives.d_multigraph.base
   :members:
   :show-inheritance:

Mixed Multi-Graph
-----------------

A mixed multi-graph contains both directed and undirected edges, and allows multiple edges
between vertices. This is the underlying graph structure for semi-directed phylogenetic
networks, where tree edges are undirected and reticulation arcs are directed.

.. automodule:: phylozoo.core.primitives.m_multigraph.base
   :members:
   :show-inheritance:
