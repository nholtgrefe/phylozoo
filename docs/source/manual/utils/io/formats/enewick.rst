eNewick
=======

Newick is the standard text format for phylogenetic trees (nested parentheses,
labels, branch lengths). **eNewick** (Extended Newick) extends it to support
directed phylogenetic networks: hybrid nodes, gamma probabilities, and edge
attributes.

PhyloZoo uses the **same eNewick format** for both directed and semi-directed
networks. The file on disk is always eNewick. The difference is how each class
handles it:

- :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` — Reads and writes eNewick directly (the network is
  already rooted).

- :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` — When **writing**, the semi-directed network is
  rooted (a root location is chosen, by default automatically), then the
  resulting directed network is written as eNewick. When **reading**, the
  eNewick string is parsed into a directed network, then converted to
  semi-directed (the root choice is discarded and undirected tree edges are
  restored). So the same eNewick file can be loaded as either class; the
  in-memory representation differs (rooted vs. unrooted).

.. seealso::
   `Newick format <https://en.wikipedia.org/wiki/Newick_format>`_ — Wikipedia (eNewick extends Newick)

Classes and extensions
----------------------

**Classes:** :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`,
:class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`

**File extensions:** ``.enewick``, ``.eNewick``, ``.enwk``, ``.nwk``, ``.newick``, ``.enw``

**Format name:** ``'enewick'`` (default for both classes)

Structure
---------

eNewick uses nested parentheses, commas, node labels, and optional branch lengths
and comments. Reticulations use ``#H1``, ``#H2``, … markers. Edge attributes such
as gamma can appear in **comments** in square brackets (e.g. ``[&gamma=0.6]``)
before the child’s branch length.

.. code-block:: text

   ((A:1.0[&gamma=0.6],B:1.0[&gamma=0.4])#H1:0.5,C:2.0);

This example shows a hybrid node ``#H1`` with two incoming edges (gamma 0.6 and 0.4)
and a branch length 0.5. Simple trees look like standard Newick, e.g. ``((A,B),C);``.

Examples
--------

**Directed network (read/write eNewick directly):**

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

**Semi-directed network (same eNewick; roots when writing, converts when reading):**

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
   network.save("network.enewick")   # rooted, then written as eNewick
   network2 = SemiDirectedPhyNetwork.load("network.enewick")  # parsed to D, then D→SD

You can pass ``root_location`` when saving a semi-directed network to control
where the root is placed before writing eNewick.

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`dot` — DOT and PhyloZoo-DOT for networks
- :doc:`../../../core/networks/semi_directed/overview` — Semi-directed networks
- :doc:`../../../core/networks/directed/overview` — Directed networks
