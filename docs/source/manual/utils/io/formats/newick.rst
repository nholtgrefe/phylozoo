Newick and eNewick Formats
==========================

Newick is the standard text format for phylogenetic trees. PhyloZoo uses **standard
Newick** for semi-directed networks (trees and graphs with undirected edges) and
**eNewick** (Extended Newick) for directed phylogenetic networks with hybrid nodes
and edge attributes.

Standard Newick
---------------

Standard Newick represents trees (and tree-like networks) with nested parentheses,
commas, labels, and optional branch lengths.

**Classes using Newick format:** :class:`phylozoo.core.network.sdnetwork.SemiDirectedPhyNetwork` (default format).

**File extensions:** ``.nwk``, ``.newick``, ``.enewick``, ``.eNewick``, ``.enw``

Structure:

.. code-block:: text

   ((A,B),C);

For semi-directed networks, undirected edges are represented in a way compatible with
standard Newick parsers.

Example:

.. code-block:: python

   from phylozoo import SemiDirectedPhyNetwork
   network = SemiDirectedPhyNetwork(
       undirected_edges=[(3, 1), (3, 2), (3, 4)],
       nodes=[
           (1, {"label": "A"}),
           (2, {"label": "B"}),
           (4, {"label": "C"})
       ]
   )
   network.save("network.newick")
   network2 = SemiDirectedPhyNetwork.load("network.newick")

eNewick (Extended Newick)
-------------------------

eNewick extends Newick to support directed networks with hybrid nodes and edge
attributes such as branch lengths and gamma probabilities.

**Classes using eNewick format:** :class:`phylozoo.core.network.dnetwork.DirectedPhyNetwork` (default format).

**File extensions:** ``.enewick``, ``.eNewick``, ``.enwk``, ``.nwk``, ``.newick``

Structure:

.. code-block:: text

   ((A:1.0[&gamma=0.6],B:1.0[&gamma=0.4])#H1:0.5,C:2.0);

This example shows a hybrid node ``#H1`` with two incoming edges (gamma 0.6 and 0.4)
and a branch length 0.5.

Example:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   network = DirectedPhyNetwork(
       edges=[
           ("root", "u1"), ("root", "u2"),
           ("u1", "h", {"gamma": 0.6}),
           ("u2", "h", {"gamma": 0.4}),
           ("h", "leaf1")
       ],
       nodes=[("leaf1", {"label": "A"})]
   )
   network.save("network.enewick")
   network2 = DirectedPhyNetwork.load("network.enewick")

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`dot` — DOT and PhyloZoo-DOT for networks
   * :doc:`../../../core/networks/semi_directed/overview` — Semi-directed networks
   * :doc:`../../../core/networks/directed/overview` — Directed networks
