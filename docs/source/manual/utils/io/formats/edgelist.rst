Edge list
=========

Edge list is a simple text format for graphs: one edge per line, given as two node
identifiers. PhyloZoo uses it for directed multigraphs (e.g. graph topology without
labels or attributes).

.. seealso::
   `Edge list <https://en.wikipedia.org/wiki/Edge_list>`_ — Wikipedia

Classes and extensions
----------------------

**Classes:** :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`

**File extensions:** ``.edgelist``, ``.edges``

Structure
---------

Each line contains an edge as two node identifiers (space- or tab-separated):

.. code-block:: text

   1 2
   2 3
   3 1

Examples
--------

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   graph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (3, 1)])
   graph.save("graph.edgelist", format="edgelist")
   graph2 = DirectedMultiGraph.load("graph.edgelist", format="edgelist")

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`dot` — DOT format for directed graphs with attributes
- :doc:`../../../core/primitives/directed_multigraph` — Directed multigraphs
