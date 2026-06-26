DOT
===

DOT is the Graphviz graph-description language. PhyloZoo uses it in two flavours:

- **DOT** (``dot``): standard Graphviz, parseable by any DOT tool.
- **PhyloZoo-DOT** (``phylozoo-dot``): a compact PhyloZoo dialect for
  semi-directed/mixed structures.

.. seealso::
   `DOT (graph description language) <https://en.wikipedia.org/wiki/DOT_(graph_description_language)>`_ — Wikipedia

Classes and extensions
----------------------

**DOT** (``.dot``, ``.gv``):

- :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` (default),
  :class:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator` (default),
  :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` — directed only.
- :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`,
  :class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator` — mixed
  graphs, with undirected edges encoded via ``dir=none`` (see below).

**PhyloZoo-DOT** (``.pzdot``):

- :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` (default),
  :class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator` (default),
  :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`.

Generators (:class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator`,
:class:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator`) inherit the
formats of their underlying graph, so they can be saved and loaded directly — no need
to route through ``generator.graph``.

Parallel (multi-)edges are encoded in both flavours with an explicit
``key=<int>`` attribute.

Standard DOT
------------

Standard DOT uses ``digraph`` and ``->``. A purely directed graph writes directly:

.. code-block:: text

   digraph {
       root -> u1;
       root -> u2;
       u1 -> h [label="gamma=0.6"];
   }

A **mixed** graph is still written as a ``digraph``; its **undirected** edges
carry ``dir=none`` — Graphviz's own marker for an edge drawn without arrowheads.
The file is therefore valid DOT (it opens in any Graphviz tool) and round-trips
losslessly, since ``dir=none`` is read back as an undirected edge:

.. code-block:: text

   digraph {
       1 -> 2;                 // directed
       1 -> 2 [key=1];         // parallel directed edge
       2 -> 3 [dir=none];      // undirected
   }

PhyloZoo-DOT
------------

PhyloZoo-DOT instead uses a ``graph`` block and distinguishes edge types
syntactically — ``->`` for directed and ``--`` for undirected — which is compact
and human-readable:

.. code-block:: text

   graph {
       1 -> 2;                 // directed
       2 -- 3;                 // undirected
       3 -- 4 [key=1];         // parallel undirected edge
   }

.. note::

   PhyloZoo-DOT mixes ``->`` and ``--`` in one ``graph`` block, which is **not**
   valid Graphviz DOT (a ``graph`` may not contain ``->``). Use the standard
   ``dot`` format if the file must be read by Graphviz or another DOT tool; use
   ``phylozoo-dot`` for compact PhyloZoo-internal storage.

Examples
--------

**Mixed multigraph, both flavours:**

.. code-block:: python

   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

   G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])

   G.save("graph.pzdot")                       # phylozoo-dot (default)
   G.save("graph.dot", format="dot")           # standard, Graphviz-readable

   G2 = MixedMultiGraph.load("graph.dot")       # auto-detects 'dot' from .dot
   G3 = MixedMultiGraph.load("graph.pzdot")     # auto-detects 'phylozoo-dot'

**Convert between the two flavours:**

.. code-block:: python

   MixedMultiGraph.convert("graph.pzdot", "graph.dot")

**Generators (saved/loaded directly):**

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import (
       SemiDirectedGenerator, all_level_k_generators,
   )

   gen = next(all_level_k_generators(3))
   gen.save("gen.pzdot")                        # phylozoo-dot (default)
   gen.save("gen.dot", format="dot")            # standard Graphviz
   gen2 = SemiDirectedGenerator.load("gen.dot")  # -> SemiDirectedGenerator

**Directed network (standard DOT):**

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B")],
       nodes=[("A", {"label": "A"}), ("B", {"label": "B"})]
   )
   network.save("network.dot", format="dot")
   network2 = DirectedPhyNetwork.load("network.dot", format="dot")

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`enewick` — eNewick for trees and networks
- :doc:`../../../core/networks/semi_directed/overview` — Semi-directed networks
- :doc:`../../../core/primitives/directed_multigraph` — Directed multigraphs
- :doc:`../../../core/primitives/mixed_multigraph` — Mixed multigraphs
- :doc:`../../../visualization/overview` — Visualization (plotting, styling)
