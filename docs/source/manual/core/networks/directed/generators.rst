Directed Network Generators
============================

The :mod:`phylozoo.core.network.dnetwork.generator` module provides the :class:`DirectedGenerator`
class and related functions for working with level-k network generators. Generators are minimal
biconnected components that represent the core structure of level-k directed phylogenetic networks
:cite:`PhyloZoo2024`. They are used to characterize and construct networks based on their level.

All classes and functions on this page can be imported from the core network generator module:

.. code-block:: python

   from phylozoo.core.dnetwork.generator import DirectedGenerator, all_level_k_generators, generators_from_network
   # or directly
   from phylozoo.core.dnetwork.generator import DirectedGenerator

Working with Directed Generators
---------------------------------

Directed generators represent the fundamental building blocks of level-k directed phylogenetic
networks. Each generator is a minimal biconnected component that can be combined with other
generators to construct networks of a specific level.

Creating Generators
^^^^^^^^^^^^^^^^^^^

Directed generators are minimal biconnected components that represent the core structure
of level-k directed phylogenetic networks. They are created from DirectedMultiGraph objects
that represent the generator topology. Each generator has a level (k), hybrid nodes, and
sides (attachment points) where networks can be connected.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import DirectedGenerator, all_level_k_generators
   from phylozoo.core.primitives.dmultigraph import DirectedMultiGraph
   
   # Create a level-1 generator
   gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])  # Parallel edges
   generator = DirectedGenerator(gen_graph)
   
   # Access properties
   level = generator.level  # 1
   hybrid_nodes = generator.hybrid_nodes  # {4}
   sides = generator.sides  # List of attachment points

Accessing Generator Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generators provide properties to access their structure and components:

.. code-block:: python

   # Access generator properties
   level = generator.level           # The level (k) of the generator
   hybrid_nodes = generator.hybrid_nodes  # Set of hybrid nodes
   sides = generator.sides           # List of sides (attachment points)

Generator Operations
--------------------

Extracting Generators from Networks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.dnetwork.generator.generators_from_network` function extracts
all generators from a binary directed network. This function decomposes the network
into its biconnected components and identifies the generator structures. The network
must be binary (no parallel edges or high-degree nodes) for accurate generator extraction.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import generators_from_network
   from phylozoo import DirectedPhyNetwork
   
   # Extract generators from a network
   network = DirectedPhyNetwork.load("network.enewick")
   generators = list(generators_from_network(network))
   
   for gen in generators:
       print(f"Level: {gen.level}, Hybrid nodes: {gen.hybrid_nodes}")

Generating All Level-k Generators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.dnetwork.generator.all_level_k_generators` function generates
all possible level-k directed generators. This is useful for network classification,
as any level-k network can be constructed by combining these generators. The function
returns a set of all distinct generator structures for the specified level.

.. code-block:: python

   # Get all level-2 generators
   level_2_generators = all_level_k_generators(2)
   print(f"Number of level-2 generators: {len(level_2_generators)}")

Generator Sides
^^^^^^^^^^^^^^^

Generators have "sides" which are attachment points where networks can be connected
to form larger networks. There are two types of sides: HybridSide (attachment at a
hybrid node) and DirEdgeSide (attachment along a directed edge). Sides define how
generators can be combined to construct level-k networks.

**HybridSide**

The :class:`phylozoo.core.dnetwork.generator.HybridSide` class represents attachment
at a hybrid node.

**DirEdgeSide**

The :class:`phylozoo.core.dnetwork.generator.DirEdgeSide` class represents attachment
along a directed edge.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import Side, HybridSide, DirEdgeSide
   
   # Sides represent attachment points
   # HybridSide: attachment at a hybrid node
   # DirEdgeSide: attachment along a directed edge
   
   generator = DirectedGenerator(...)
   for side in generator.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, DirEdgeSide):
           print(f"Edge side: {side.edge}")

.. note::
   Generators are used internally for network classification and construction. They 
   are particularly useful for understanding the structure of level-k networks.

.. tip::
   Use ``all_level_k_generators(k)`` to enumerate all possible level-k generator 
   structures. This is useful for network classification and generation.

.. warning::
   Extracting generators from networks requires the network to be binary. Networks 
   with parallel edges may produce unexpected results.

See Also
--------

- :doc:`Directed Networks (Advanced) <advanced>` - Network level classification
- :doc:`Semi-Directed Network Generators <../semi_directed/generators>` - Semi-directed generators
