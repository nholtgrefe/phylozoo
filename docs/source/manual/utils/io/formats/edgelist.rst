Edge List Format
================

Edge list is a simple text format for graphs: one edge per line, given as two node
identifiers. PhyloZoo uses it for directed multigraphs (e.g. graph topology without
labels or attributes).

**Classes using Edge list format:** :class:`phylozoo.core.primitives.d_multigraph.DirectedMultiGraph`

**File extensions:** ``.edgelist``, ``.edges``

Format Structure
----------------

Each line contains an edge as two node identifiers (space- or tab-separated):

.. code-block:: text

   1 2
   2 3
   3 1

Example
-------

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   graph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (3, 1)])
   graph.save("graph.edgelist", format="edgelist")
   graph2 = DirectedMultiGraph.load("graph.edgelist", format="edgelist")

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`dot` — DOT format for directed graphs with attributes
   * :doc:`../../../core/primitives/overview` — Primitives (multigraphs)
