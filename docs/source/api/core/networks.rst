Networks
========

The network module provides classes for representing phylogenetic networks. PhyloZoo supports
two main types of networks: directed phylogenetic networks (rooted networks with directed edges)
and semi-directed phylogenetic networks (networks with both directed and undirected edges).

Directed Phylogenetic Networks
-------------------------------

Directed phylogenetic networks are rooted networks where all edges are directed. Internal nodes
are either tree nodes (in-degree 1, out-degree :math:`\geq 2`) or hybrid nodes (in-degree
:math:`\geq 2`, out-degree 1).

.. automodule:: phylozoo.core.network.dnetwork
   :members:
   :show-inheritance:

Semi-Directed Phylogenetic Networks
-----------------------------------

Semi-directed phylogenetic networks contain both directed edges (reticulation arcs) and
undirected edges (tree edges). Internal nodes have degree :math:`\geq 3`, and hybrid nodes
have at least one incoming directed edge.

.. automodule:: phylozoo.core.network.sdnetwork
   :members:
   :show-inheritance:
