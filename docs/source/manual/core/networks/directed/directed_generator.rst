Directed Generator
==================

The :mod:`phylozoo.core.network.dnetwork.generator` module provides the :class:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator`
class and related functions for working with directed level-k generators (sometimes also referred to as "rooted level-k generators").

All classes and functions on this page can be imported from the directed generator module:

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import *
   # or directly
   from phylozoo.core.network.dnetwork.generator import DirectedGenerator

What is a (directed) generator?
----------------------------------------

A **directed generator** is a fundamental building block of a binary directed phylogenetic network.
Formally, it is a directed multi-graph that can be obtained from a blob of a directed network 
by suppressing all vertices with degree 2 except those with indegree 2.
Equivalently, a directed generator can be obtained from a simple directed 
network (i.e., a network consisting of a single internal blob) by removing all 
leaves and thereafter suppressing all vertices with degree 2 except those with indegree 2.
See, for example, :cite:`Gambette2009` for more details.

Working with generators
-----------------------

Creating a generator
^^^^^^^^^^^^^^^^^^^^

You can obtain a directed generator in two ways.

**From a graph**

Build a generator from a :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`
that represents the generator topology. A level-0 generator is a single node with no edges. A level-1 generator is a root
node with two parallel directed edges to one hybrid node.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import DirectedGenerator
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   # Level-1 generator (root with two parallel edges to hybrid node)
   gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
   generator = DirectedGenerator(gen_graph)
   
   # Level-0 generator (single node, no edges)
   gen_graph0 = DirectedMultiGraph()
   gen_graph0.add_node(1)
   generator0 = DirectedGenerator(gen_graph0)

**From a network**

Alternatively, you can extract all generators from a directed network using the :func:`~phylozoo.core.network.dnetwork.generator.base.generators_from_network` function.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import generators_from_network
   from phylozoo import DirectedPhyNetwork
   
   network = DirectedPhyNetwork.load("network.enewick")
   generators = list(generators_from_network(network))
   for gen in generators:
       print(f"Level: {gen.level}, Hybrid nodes: {gen.hybrid_nodes}")

.. note::
   Upon construction, a generator is validated using :meth:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator.validate` to ensure it has a valid structure that adheres to the definition of a directed generator.
   For details on PhyloZoo's validation system and how to disable it for performance reasons, see
   :doc:`Validation <../../../utils/validation>`.

Accessing generator attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The underlying graph, root node, level, and hybrid nodes of a generator can be accessed as follows.
The hybrid nodes are the nodes with indegree 2; the level is the number of hybrid nodes; the root
node is the node with in-degree 0.

.. code-block:: python

   graph = generator.graph
   root = generator.root_node
   level = generator.level
   hybrid_nodes = generator.hybrid_nodes

Sides of a generator
^^^^^^^^^^^^^^^^^^^^

The sides of a generator are precisely the points where leaves can be attached to create a directed network.
The module provides the following types of side types.

**Base classes**

- :class:`~phylozoo.core.network.dnetwork.generator.side.Side` — The base class from which all other side types inherit.
- :class:`~phylozoo.core.network.dnetwork.generator.side.NodeSide` is the base class for node attachment.
- :class:`~phylozoo.core.network.dnetwork.generator.side.EdgeSide` is the base class for edge attachment.

**Node sides**

- :class:`~phylozoo.core.network.dnetwork.generator.side.HybridSide` — A hybrid node (in-degree ≥ 2) with out-degree 0 that serves as an attachment point.
- :class:`~phylozoo.core.network.dnetwork.generator.side.IsolatedNodeSide` — The single vertex of a level-0 generator. Exactly one such side exists when level is 0.

**Edge sides**

- :class:`~phylozoo.core.network.dnetwork.generator.side.DirEdgeSide` — A directed edge (possibly one of several parallel edges), identified by ``u``, ``v``, and ``key``.

The sides of a generator can be accessed as follows:

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import (
       HybridSide, DirEdgeSide, IsolatedNodeSide,
   )
   
   hybrid_side = HybridSide(node=3)
   edge_side = DirEdgeSide(u=8, v=4, key=0)
   
   # generator.sides: level-0 → [IsolatedNodeSide(root)]; level≥1 → edge_sides + hybrid_sides
   all_sides = generator.sides
   hybrid_sides = generator.hybrid_sides
   edge_sides = generator.edge_sides
   parallel_sides = generator.parallel_edge_sides
   non_parallel_sides = generator.non_parallel_edge_sides
   
   for side in generator.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, DirEdgeSide):
           print(f"Edge side: ({side.u}, {side.v}, key={side.key})")
       elif isinstance(side, IsolatedNodeSide):
           print(f"Isolated node side: {side.node}")

Attaching leaves to a generator
--------------------------------

:func:`~phylozoo.core.network.dnetwork.generator.attachment.attach_leaves_to_generator` takes a
generator and a mapping ``side_taxa`` from sides to ordered lists of taxon labels, and returns a
:class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` by attaching new leaf nodes to
the generator's graph.

**Rules**

- Every :class:`~phylozoo.core.network.dnetwork.generator.side.HybridSide` must appear in
  ``side_taxa`` with **exactly one** taxon.
- Every :class:`~phylozoo.core.network.dnetwork.generator.side.IsolatedNodeSide` (level-0 only)
  must appear in ``side_taxa`` with **exactly three** taxa.
- :class:`~phylozoo.core.network.dnetwork.generator.side.DirEdgeSide` sides may be omitted or
  given an empty list; they then receive no leaves. Note that the order of the taxa on the sides is important.
- At least two taxa must be attached in total across all sides.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import (
       DirectedGenerator, attach_leaves_to_generator,
   )
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   dmg = DirectedMultiGraph(edges=[(0, 1), (0, 1)])
   gen = DirectedGenerator(dmg)
   hybrid_side = gen.hybrid_sides[0]
   edge_side = gen.edge_sides[0]
   network = attach_leaves_to_generator(
       gen, {hybrid_side: ["H"], edge_side: ["A", "B"]}
   )
   sorted(network.taxa)
   # ['A', 'B', 'H']

Enumerating all generators
--------------------------

The function :func:`~phylozoo.core.network.dnetwork.generator.construction.all_level_k_generators` generates
all distinct level-k directed generators (up to isomorphism), using the algorithm described in :cite:`Gambette2009`. It starts with level-0 (single node),
then iteratively applies two extension rules (R1 and R2) to level-(k-1) generators and
deduplicates by isomorphism.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import all_level_k_generators
   
   level_0 = all_level_k_generators(0)   # 1 generator
   level_1 = all_level_k_generators(1)    # 1 generator
   level_2 = all_level_k_generators(2)    # 4 generators

See Also
--------

- :doc:`Directed Network Class <directed_network_class>` - The DirectedPhyNetwork class
- :doc:`Semi-Directed Network Generators <../semi_directed/generators>` - Semi-directed generators
- :doc:`Validation <../../../utils/validation>` - Validation system and disabling validation
- :doc:`API Reference <../../../../api/core/network/index>` - Complete function signatures and detailed examples
