Directed Generator
==================

The :mod:`phylozoo.core.network.dnetwork.generator` module provides the :class:`DirectedGenerator`
class and related functions for working with level-k network generators. Generators are minimal
biconnected components that represent the core structure of level-k directed phylogenetic networks
:cite:`Gambette2009`. They are used to characterize and construct networks based on their level.

All classes and functions on this page can be imported from the core network generator module:

.. code-block:: python

   from phylozoo.core.dnetwork.generator import (
       DirectedGenerator,
       all_level_k_generators,
       generators_from_network,
       Side,
       HybridSide,
       DirEdgeSide,
   )

DirectedGenerator Class
-----------------------

The :class:`phylozoo.core.dnetwork.generator.DirectedGenerator` class represents a level-k generator
for directed phylogenetic networks. A generator is a biconnected component that represents the
core structure of a level-k network. Generators use DirectedMultiGraph directly and have their
own validation rules. These are not networks themselves, but simplified structures used to build networks.

Creating Generators
^^^^^^^^^^^^^^^^^^^

Generators are created from :class:`phylozoo.core.primitives.d_multigraph.DirectedMultiGraph` objects
that represent the generator topology. Each generator has a level (k), hybrid nodes, and
sides (attachment points) where networks can be connected.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import DirectedGenerator
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   # Create a level-1 generator (parallel edges from root to hybrid node)
   gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])  # Parallel edges
   generator = DirectedGenerator(gen_graph)
   
   # Create a level-0 generator (single node)
   gen_graph0 = DirectedMultiGraph()
   gen_graph0.add_node(1)
   generator0 = DirectedGenerator(gen_graph0)

Accessing Generator Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`DirectedGenerator` class provides read-only access to fundamental structural properties:

**Basic Properties**

.. code-block:: python

   # Get the underlying graph structure
   graph = generator.graph  # Returns DirectedMultiGraph
   
   # Get the root node (unique node with in-degree 0)
   root = generator.root_node  # Returns the root node ID

**Level and Hybrid Nodes**

.. code-block:: python

   # Get the level (number of hybrid nodes)
   level = generator.level  # Returns int, e.g., 1
   
   # Get all hybrid nodes (nodes with in-degree >= 2)
   hybrid_nodes = generator.hybrid_nodes  # Returns set of node IDs

**Generator Sides**

Sides are attachment points where networks can be connected to form larger networks. There are two
types of sides: :class:`HybridSide` (attachment at a hybrid node) and :class:`DirEdgeSide`
(attachment along a directed edge).

.. code-block:: python

   # Get all sides (both edge sides and hybrid sides)
   all_sides = generator.sides  # Returns list[Side]
   
   # Get only hybrid sides (hybrid nodes with out-degree 0)
   hybrid_sides = generator.hybrid_sides  # Returns list[HybridSide]
   
   # Get all edge sides (all edges as DirEdgeSide objects)
   edge_sides = generator.edge_sides  # Returns list[DirEdgeSide]
   
   # Get parallel edge sides (groups of parallel edges)
   parallel_sides = generator.parallel_edge_sides  # Returns list[tuple[DirEdgeSide, ...]]
   
   # Get non-parallel edge sides
   non_parallel_sides = generator.non_parallel_edge_sides  # Returns list[DirEdgeSide]

**Validation**

Generators are validated at construction time. The :meth:`validate` method can be called
explicitly to check the generator structure:

.. code-block:: python

   # Validate generator structure
   generator.validate()  # Raises exception if generator is invalid

Generator Functions
-------------------

The generator module provides functions to extract generators from networks and generate all
possible level-k generators.

**Extracting Generators from Networks**

The :func:`phylozoo.core.dnetwork.generator.generators_from_network` function extracts
all generators from a binary directed network. This function decomposes the network
into its biconnected components and identifies the generator structures.

The network must be binary (no parallel edges or high-degree nodes) for accurate generator extraction.
If the network has parallel edges, a warning is issued but the function continues.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import generators_from_network
   from phylozoo import DirectedPhyNetwork
   
   # Extract generators from a network
   network = DirectedPhyNetwork.load("network.enewick")
   generators = list(generators_from_network(network))
   
   for gen in generators:
       print(f"Level: {gen.level}, Hybrid nodes: {gen.hybrid_nodes}")

**Generating All Level-k Generators**

The :func:`phylozoo.core.dnetwork.generator.all_level_k_generators` function generates
all possible level-k directed generators. This is useful for network classification,
as any level-k network can be constructed by combining these generators. The function
returns a set of all distinct generator structures for the specified level.

The function constructs generators by starting with level-0 generators and iteratively
applying R1 and R2 rules. Validation is deferred during construction for performance
optimization, but generators can be explicitly validated using the :meth:`validate` method.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import all_level_k_generators
   
   # Get all level-0 generators (single node)
   level_0_generators = all_level_k_generators(0)
   print(f"Number of level-0 generators: {len(level_0_generators)}")  # 1
   
   # Get all level-1 generators (parallel edges)
   level_1_generators = all_level_k_generators(1)
   print(f"Number of level-1 generators: {len(level_1_generators)}")  # 1
   
   # Get all level-2 generators
   level_2_generators = all_level_k_generators(2)
   print(f"Number of level-2 generators: {len(level_2_generators)}")  # 4

Generator Sides
---------------

Generator sides represent attachment points where additional structure (leaves, trees) can be
attached when building generators from lower-level generators. The module provides three classes
for representing sides.

**Side Base Class**

The :class:`phylozoo.core.dnetwork.generator.Side` class is the base class for all sides.
It serves as a general wrapper class. Use :class:`HybridSide` for node sides and
:class:`DirEdgeSide` for edge sides.

**HybridSide**

The :class:`phylozoo.core.dnetwork.generator.HybridSide` class represents a hybrid side of a generator.
A hybrid side is a node with in-degree 2 and out-degree 0, representing a hybrid node that serves
as an attachment point.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import HybridSide
   
   # Create a hybrid side
   hybrid_side = HybridSide(node=3)
   print(hybrid_side.node)  # 3

**DirEdgeSide**

The :class:`phylozoo.core.dnetwork.generator.DirEdgeSide` class represents a directed edge side of a generator.
A directed edge side represents an edge (possibly parallel) that serves as an attachment point.
The edge is identified by its endpoints and key.

.. code-block:: python

   from phylozoo.core.dnetwork.generator import DirEdgeSide
   
   # Create an edge side
   edge_side = DirEdgeSide(u=8, v=4, key=0)
   print(edge_side.u)    # 8
   print(edge_side.v)    # 4
   print(edge_side.key)  # 0

**Working with Sides**

You can iterate over generator sides and check their types:

.. code-block:: python

   from phylozoo.core.dnetwork.generator import HybridSide, DirEdgeSide
   
   generator = DirectedGenerator(...)
   for side in generator.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, DirEdgeSide):
           print(f"Edge side: ({side.u}, {side.v}, key={side.key})")

.. note::
   Generators are used internally for network classification and construction. They 
   are particularly useful for understanding the structure of level-k networks.

.. tip::
   Use ``all_level_k_generators(k)`` to enumerate all possible level-k generator 
   structures. This is useful for network classification and generation.

.. warning::
   Extracting generators from networks requires the network to be binary. Networks 
   with parallel edges may produce unexpected results, and a warning will be issued.

See Also
--------

- :doc:`Directed Network Algorithms <directed_network_algorithms>` - Network level classification
- :doc:`Directed Network Class <directed_network_class>` - The DirectedPhyNetwork class
- :doc:`Semi-Directed Network Generators <../semi_directed/semi_directed_generator>` - Semi-directed generators
- :doc:`API Reference <../../../api/core/network>` - Complete function signatures and detailed examples