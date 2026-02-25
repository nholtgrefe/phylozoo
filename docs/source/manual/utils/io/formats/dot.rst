DOT and PhyloZoo-DOT Formats
============================

DOT is the Graphviz format for graph visualization. PhyloZoo supports **standard DOT**
for directed graphs and **PhyloZoo-DOT** (``pzdot``) for semi-directed networks,
where the distinction between directed and undirected edges is preserved.

DOT Format
----------

Standard DOT uses a declarative syntax for directed graphs (``digraph``).

**Classes using DOT format:**

* :class:`phylozoo.core.network.dnetwork.DirectedPhyNetwork`
* :class:`phylozoo.core.primitives.d_multigraph.DirectedMultiGraph` (default format)

**File extensions:** ``.dot``, ``.gv``

Structure:

.. code-block:: text

   digraph {
       root -> u1;
       root -> u2;
       u1 -> h [label="gamma=0.6"];
       u2 -> h [label="gamma=0.4"];
       h -> leaf1 [label="A"];
   }

Example:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B")],
       nodes=[("A", {"label": "A"}), ("B", {"label": "B"})]
   )
   network.save("network.dot", format="dot")
   network2 = DirectedPhyNetwork.load("network.dot", format="dot")

PhyloZoo-DOT Format
-------------------

PhyloZoo-DOT is a custom variant for semi-directed networks. It uses ``graph`` (not
``digraph``) and represents undirected edges with ``--`` and directed edges with
``->``, so that both edge types can be represented in one file.

**Classes using PhyloZoo-DOT format:**

* :class:`phylozoo.core.network.sdnetwork.SemiDirectedPhyNetwork`
* :class:`phylozoo.core.primitives.m_multigraph.MixedMultiGraph` (default format)

**File extensions:** ``.pzdot``

Structure:

.. code-block:: text

   graph {
       node1 -- node2 [dir=none];
       node2 -> node3;
   }

Example:

.. code-block:: python

   from phylozoo import SemiDirectedPhyNetwork
   network = SemiDirectedPhyNetwork(
       directed_edges=[(5, 4, {"gamma": 0.6})],
       undirected_edges=[(4, 1), (4, 2)],
       nodes=[
           (1, {"label": "A"}),
           (2, {"label": "B"})
       ]
   )
   network.save("network.pzdot", format="phylozoo-dot")
   network2 = SemiDirectedPhyNetwork.load("network.pzdot", format="phylozoo-dot")

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`newick` — Newick and eNewick for trees and networks
   * :doc:`../../../core/networks/semi_directed/overview` — Semi-directed networks
   * :doc:`../../../visualization/overview` — Visualization (layouts, plotting)
