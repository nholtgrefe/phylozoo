Overview
========

This section covers working with phylogenetic networks in PhyloZoo. PhyloZoo provides 
two main network classes for representing evolutionary relationships.

Network Types
-------------

**DirectedPhyNetwork**
   Fully directed phylogenetic networks with a single root node. All edges are directed, 
   explicitly representing the direction of evolutionary time. Hybrid nodes have 
   in-degree >= 2 and out-degree 1.
   
   **I/O Formats**: eNewick (default), DOT. See :doc:`I/O <../../io>` for details.
   
   See :doc:`Directed Networks <directed/overview>` for detailed documentation.

**SemiDirectedPhyNetwork**
   Networks with directed hybrid edges and undirected tree edges. These allow for more 
   flexible representation and are useful for unrooted analyses. Semi-directed networks 
   cannot have undirected cycles.
   
   **I/O Formats**: Newick (default), PhyloZoo-DOT. See :doc:`I/O <../../io>` for details.
   
   See :doc:`Semi-Directed Networks <semi_directed/overview>` for detailed documentation.

Choosing a Network Type
-----------------------

* Use **DirectedPhyNetwork** when you need explicit root information and want to represent 
  time explicitly through directed edges.
* Use **SemiDirectedPhyNetwork** when working with unrooted analyses or when the root 
  location is uncertain.

Both network types support:
* Hybrid nodes (reticulations)
* Branch lengths
* Node and edge attributes
* Network analysis and transformation operations

For detailed information, see:

* :doc:`Directed Networks <directed/overview>`: Creating, manipulating, and analyzing directed networks
* :doc:`Semi-Directed Networks <semi_directed/overview>`: Creating, manipulating, and analyzing semi-directed networks

.. seealso::
   For I/O operations, see :doc:`I/O <../../io>`. 
   For visualization, see :doc:`Visualization <../../visualization/viz>`.
