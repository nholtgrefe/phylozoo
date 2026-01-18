Overview
========

**SemiDirectedPhyNetwork** represents networks with directed hybrid edges and undirected 
tree edges. These allow for more flexible representation and are useful for unrooted 
analyses. Semi-directed networks cannot have undirected cycles.

Key Features
------------

* **Flexible Rooting**: No fixed root; can be rooted at various locations
* **Mixed Edges**: Directed hybrid edges and undirected tree edges
* **Hybrid Nodes**: Support for reticulation events
* **Edge Attributes**: Branch lengths and other attributes
* **I/O Formats**: Newick (default), PhyloZoo-DOT

I/O Formats
-----------

* **Newick** (default): Standard Newick format. Extensions: ``.nwk``, ``.newick``, ``.enewick``, ``.eNewick``, ``.enw``
* **PhyloZoo-DOT**: Custom DOT format. Extension: ``.pzdot``

See :doc:`I/O <../../../../io>` for details.

Documentation
-------------

* :doc:`Basic Operations <basic>`: Creating networks, basic properties, loading/saving
* :doc:`Advanced Features <advanced>`: Network features, transformations, classifications
* :doc:`Generators <generators>`: Level-k network generators

.. seealso::
   For directed networks, see :doc:`Directed Networks <../directed/overview>`. 
   For I/O operations, see :doc:`I/O <../../../../io>`.
