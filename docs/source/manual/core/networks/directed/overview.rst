Overview
========

**DirectedPhyNetwork** represents fully directed phylogenetic networks with a single root 
node. All edges are directed, explicitly representing the direction of evolutionary time. 
Hybrid nodes have in-degree >= 2 and out-degree 1.

Key Features
------------

* **Single Root**: Every network has exactly one root node
* **Fully Directed**: All edges are directed, representing time flow
* **Hybrid Nodes**: Support for reticulation events (in-degree >= 2, out-degree 1)
* **Edge Attributes**: Branch lengths, bootstrap values, gamma probabilities
* **I/O Formats**: eNewick (default), DOT

I/O Formats
-----------

* **eNewick** (default): Extended Newick format. Extensions: ``.enewick``, ``.eNewick``, ``.enwk``, ``.nwk``, ``.newick``
* **DOT**: Graphviz format. Extensions: ``.dot``, ``.gv``

See :doc:`I/O <../../../../io>` for details.

Documentation
-------------

* :doc:`Basic Operations <basic>`: Creating networks, basic properties, loading/saving
* :doc:`Advanced Features <advanced>`: Network features, transformations, classifications, isomorphism
* :doc:`Generators <generators>`: Level-k network generators

.. seealso::
   For semi-directed networks, see :doc:`Semi-Directed Networks <../semi_directed/overview>`. 
   For I/O operations, see :doc:`I/O <../../../../io>`.
