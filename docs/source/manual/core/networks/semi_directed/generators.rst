Semi-Directed Generator
========================

The :mod:`phylozoo.core.network.sdnetwork.generator` module provides the
:class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator` class and related
functions for working with level-k semi-directed network generators.

Most classes and functions on this page can be imported from the semi-directed generator module:

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import *
   # or directly
   from phylozoo.core.network.sdnetwork.generator import SemiDirectedGenerator

Some side classes are defined in the directed generator module and reused here.


What is a (semi-directed) generator?
------------------------------------

A **semi-directed generator** is a fundamental building block of a level-k semi-directed
phylogenetic network. Formally, it is a mixed multi-graph that can be obtained from a blob of a
semi-directed network by suppressing all vertices with degree 2 except those with indegree 2.
Equivalently, a semi-directed generator can be obtained from a simple semi-directed network
(i.e., a network consisting of a single internal blob) by removing all leaves and thereafter
suppressing all vertices with degree 2 except those with indegree 2. Generators are built from a
:class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`.

Working with generators
-----------------------

Creating a generator
^^^^^^^^^^^^^^^^^^^^

You can obtain a semi-directed generator in two ways.

**From a graph**

Build a generator from a :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`.
A level-0 generator is a single node with no edges. A level-1 generator is a single node with one
bidirected self-loop. Since MixedMultiGraphs do not support bidirected edges, the generator is built from a single node with one undirected self-loop.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import SemiDirectedGenerator
   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   # Level-1 from MixedMultiGraph (one node, one undirected self-loop)
   gen_graph = MixedMultiGraph(undirected_edges=[(4, 4)])
   generator = SemiDirectedGenerator(gen_graph)


**From a directed generator**

The function :func:`~phylozoo.core.network.sdnetwork.generator.construction.dgenerator_to_sdgenerator` converts
a :class:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator` to a
:class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator` by: (1) turning
all non-hybrid edges into undirected edges and keeping edges into hybrid nodes directed; (2)
suppressing the degree-2 root node. Level-0 stays one node with no edges. Level-1 becomes one node
with one undirected self-loop (bidirected edge).

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import dgenerator_to_sdgenerator
   from phylozoo.core.network.dnetwork.generator import DirectedGenerator
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   d_gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
   d_gen = DirectedGenerator(d_gen_graph)
   sd_gen = dgenerator_to_sdgenerator(d_gen)
   sd_gen.level
   # 1
   sd_gen.graph.number_of_nodes()
   # 1
   list(sd_gen.graph.undirected_edges_iter(keys=True))
   # [(4, 4, 0)]

.. note::
   Upon construction, a generator is validated using :meth:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator.validate` to ensure it has a valid structure that adheres to the definition of a semi-directed generator.
   For details on PhyloZoo's validation system and how to disable it for performance reasons, see
   :doc:`Validation <../../../utils/validation>`.

Accessing generator attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The underlying graph, level, and hybrid nodes of a semi-directed generator can be accessed as follows.
The hybrid nodes are the nodes with indegree ≥ 2 from directed edges; the level is the number of
hybrid nodes. (Semi-directed generators have no single root node.)

.. code-block:: python

   graph = generator.graph
   level = generator.level
   hybrid_nodes = generator.hybrid_nodes

Sides of a generator
^^^^^^^^^^^^^^^^^^^^

The sides of a generator are precisely the points where leaves can be attached to create a semi-directed network.
The module provides node sides (from the directed generator module) and both directed and undirected edge sides.

**Base classes**

- :class:`~phylozoo.core.network.dnetwork.generator.side.Side` — The base class from which all other side types inherit.
- :class:`~phylozoo.core.network.dnetwork.generator.side.NodeSide` is the base class for node attachment.
- :class:`~phylozoo.core.network.dnetwork.generator.side.EdgeSide` is the base class for edge attachment.

**Node sides**

- :class:`~phylozoo.core.network.dnetwork.generator.side.HybridSide` — A hybrid node (in-degree ≥ 2 from directed edges) with out-degree 0 that serves as an attachment point.
- :class:`~phylozoo.core.network.dnetwork.generator.side.IsolatedNodeSide` — The single vertex of a level-0 generator. Exactly one such side exists when level is 0.

**Edge sides**

- :class:`~phylozoo.core.network.dnetwork.generator.side.DirEdgeSide` — A directed edge (identified by ``u``, ``v``, ``key``).
- :class:`~phylozoo.core.network.sdnetwork.generator.side.UndirEdgeSide` — An undirected edge (``u``, ``v``, ``key``).
- :class:`~phylozoo.core.network.sdnetwork.generator.side.BidirectedEdgeSide` — The undirected self-loop of a level-1 generator (single node, one edge). Identified by ``node`` and ``key``.

The sides of a generator can be accessed as follows:

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import (
       SemiDirectedGenerator, UndirEdgeSide, BidirectedEdgeSide,
   )
   from phylozoo.core.network.dnetwork.generator import HybridSide, DirEdgeSide, IsolatedNodeSide
   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   mm = MixedMultiGraph(undirected_edges=[(1, 1)])
   gen = SemiDirectedGenerator(mm)
   hybrid_sides = gen.hybrid_sides
   edge_sides = gen.edge_sides          # DirEdgeSide | UndirEdgeSide | BidirectedEdgeSide
   directed_edge_sides = gen.directed_edge_sides
   undirected_edge_sides = gen.undirected_edge_sides
   parallel_directed_sides = gen.parallel_directed_edge_sides
   non_parallel_directed_sides = gen.non_parallel_directed_edge_sides
   
   for side in gen.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, DirEdgeSide):
           print(f"Directed edge side: ({side.u}, {side.v}, key={side.key})")
       elif isinstance(side, UndirEdgeSide):
           print(f"Undirected edge side: ({side.u}, {side.v}, key={side.key})")
       elif isinstance(side, BidirectedEdgeSide):
           print(f"Bidirected edge side: node={side.node}, key={side.key}")
       elif isinstance(side, IsolatedNodeSide):
           print(f"Isolated node side: {side.node}")

Attaching leaves to a generator
---------------------------------

The function :func:`~phylozoo.core.network.sdnetwork.generator.attachment.attach_leaves_to_generator` takes a
semi-directed generator and a mapping ``side_taxa`` from sides to ordered lists of taxon labels, and returns a
:class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` by attaching new leaf nodes to
the generator's graph. 

**Rules**

- Every :class:`~phylozoo.core.network.dnetwork.generator.side.HybridSide` must appear in
  ``side_taxa`` with **exactly one** taxon.
- Every :class:`~phylozoo.core.network.dnetwork.generator.side.IsolatedNodeSide` (level-0 only)
  must appear in ``side_taxa`` with **exactly three** taxa (binary networks).
- Edge sides (:class:`~phylozoo.core.network.dnetwork.generator.side.DirEdgeSide`,
  :class:`~phylozoo.core.network.sdnetwork.generator.side.UndirEdgeSide`,
  :class:`~phylozoo.core.network.sdnetwork.generator.side.BidirectedEdgeSide`) may be omitted or
  given an empty list **except** for level-1: the single
  :class:`~phylozoo.core.network.sdnetwork.generator.side.BidirectedEdgeSide` must receive
  **at least one** taxon (so the self-loop is always subdivided away). Note that the order of the taxa on the sides is important.
- At least two taxa must be attached in total across all sides.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import (
       SemiDirectedGenerator, attach_leaves_to_generator,
   )
   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   
   mm = MixedMultiGraph(undirected_edges=[(1, 1)])
   gen = SemiDirectedGenerator(mm)
   hybrid_side = gen.hybrid_sides[0]
   bidir_side = gen.edge_sides[0]
   network = attach_leaves_to_generator(
       gen, {hybrid_side: ["H"], bidir_side: ["P", "Q"]}
   )
   sorted(network.taxa)
   # ['H', 'P', 'Q']

Enumerating all generators
---------------------------

The function :func:`~phylozoo.core.network.sdnetwork.generator.construction.all_level_k_generators` generates
all distinct level-k semi-directed generators (up to isomorphism). It obtains all level-k directed
generators, converts each with :func:`~phylozoo.core.network.sdnetwork.generator.construction.dgenerator_to_sdgenerator`,
then deduplicates by isomorphism.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import all_level_k_generators
   
   level_0 = all_level_k_generators(0)   # 1 generator
   level_1 = all_level_k_generators(1)   # 1 generator
   level_2 = all_level_k_generators(2)   # 2 generators

See Also
--------

- :doc:`Semi-Directed Network Class <semi_directed_network_class>` - The SemiDirectedPhyNetwork class
- :doc:`Directed Generator <../directed/directed_generator>` - Directed generators
- :doc:`Validation <../../../utils/validation>` - Validation system and disabling validation
- :doc:`API Reference <../../../../api/core/network/index>` - Complete function signatures and detailed examples
