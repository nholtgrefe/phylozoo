DOT
===

DOT is the Graphviz format for graph visualization. PhyloZoo supports **standard DOT**
for directed graphs and **PhyloZoo-DOT** (``pzdot``) for semi-directed networks,
where the distinction between directed and undirected edges is preserved.

.. seealso::
   `DOT (graph description language) <https://en.wikipedia.org/wiki/DOT_(graph_description_language)>`_ — Wikipedia

Classes and extensions
----------------------

**DOT:** :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`,
:class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` (default format).
Extensions: ``.dot``, ``.gv``

**PhyloZoo-DOT:** :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`,
:class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` (default format).
Extensions: ``.pzdot``

Structure
---------

Both formats use a declarative syntax with node and edge declarations. Nodes can
have attributes (e.g. labels); edges can have attributes (e.g. gamma, branch length).
The main difference is that standard DOT uses ``digraph`` and ``->`` for directed
edges only, while PhyloZoo-DOT uses ``graph`` with both ``->`` (directed) and ``--``
(undirected) to represent semi-directed networks.

DOT (standard)
^^^^^^^^^^^^^^

Standard DOT uses ``digraph`` for directed graphs. Node names in the output use
labels where available for readability.

.. code-block:: text

   digraph {
       root -> u1;
       root -> u2;
       u1 -> h [label="gamma=0.6"];
       u2 -> h [label="gamma=0.4"];
       h -> leaf1 [label="A"];
   }

PhyloZoo-DOT
^^^^^^^^^^^^

PhyloZoo-DOT uses ``graph`` (not ``digraph``) and represents undirected edges with
``--`` and directed edges with ``->``, so that both edge types can be represented
in one file.

.. code-block:: text

   graph {
       node1 -- node2 [dir=none];
       node2 -> node3;
   }

Examples
--------

**DOT (directed network):**

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B")],
       nodes=[("A", {"label": "A"}), ("B", {"label": "B"})]
   )
   network.save("network.dot", format="dot")
   network2 = DirectedPhyNetwork.load("network.dot", format="dot")

**PhyloZoo-DOT (semi-directed network):**

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

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`enewick` — eNewick for trees and networks
- :doc:`../../../core/networks/semi_directed/overview` — Semi-directed networks
- :doc:`../../../core/primitives/directed_multigraph` — Directed multigraphs
- :doc:`../../../core/primitives/mixed_multigraph` — Mixed multigraphs
- :doc:`../../../visualization/overview` — Visualization (plotting, styling)
