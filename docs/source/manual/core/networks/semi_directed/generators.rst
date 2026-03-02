Semi-Directed Generator
========================

The :mod:`phylozoo.core.network.sdnetwork.generator` module provides the :class:`SemiDirectedGenerator`
class and related functions for working with level-k network generators. Generators are minimal
biconnected components that represent the core structure of level-k semi-directed phylogenetic networks
:cite:`Gambette2009`. They are used to characterize and construct networks based on their level.

All classes and functions on this page can be imported from the core network generator module:

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import (
       SemiDirectedGenerator,
       all_level_k_generators,
       dgenerator_to_sdgenerator,
       Side,
       HybridSide,
       DirEdgeSide,
       UndirEdgeSide,
   )

SemiDirectedGenerator Class
----------------------------

The :class:`phylozoo.core.network.sdnetwork.generator.SemiDirectedGenerator` class represents a level-k generator
for semi-directed phylogenetic networks. A generator is a biconnected component that represents the
core structure of a level-k network. Generators use MixedMultiGraph directly and have their
own validation rules. These are not networks themselves, but simplified structures used to build networks.

Creating Generators
^^^^^^^^^^^^^^^^^^^

Generators are created from :class:`phylozoo.core.primitives.m_multigraph.MixedMultiGraph` objects
that represent the generator topology. Each generator has a level (k), hybrid nodes, and
sides (attachment points) where networks can be connected. Semi-directed generators can also be
created by converting directed generators using :func:`dgenerator_to_sdgenerator`.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import SemiDirectedGenerator, dgenerator_to_sdgenerator
   from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
   from phylozoo.core.network.dnetwork.generator import DirectedGenerator
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   # Create a level-1 generator (parallel directed edges from root to hybrid node)
   gen_graph = MixedMultiGraph(
       directed_edges=[(8, 4), (8, 4)],  # Parallel directed edges
       undirected_edges=[]
   )
   generator = SemiDirectedGenerator(gen_graph)
   
   # Or convert from a directed generator
   d_gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
   d_gen = DirectedGenerator(d_gen_graph)
   sd_gen = dgenerator_to_sdgenerator(d_gen)

Accessing Generator Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`SemiDirectedGenerator` class provides read-only access to fundamental structural properties:

**Basic Properties**

.. code-block:: python

   # Get the underlying graph structure
   graph = generator.graph  # Returns MixedMultiGraph

**Level and Hybrid Nodes**

.. code-block:: python

   # Get the level (number of hybrid nodes)
   level = generator.level  # Returns int, e.g., 1
   
   # Get all hybrid nodes (nodes with in-degree >= 2 from directed edges)
   hybrid_nodes = generator.hybrid_nodes  # Returns set of node IDs

**Generator Sides**

Sides are attachment points where networks can be connected to form larger networks. There are three
types of sides: :class:`HybridSide` (attachment at a hybrid node), :class:`DirEdgeSide`
(attachment along a directed edge), and :class:`UndirEdgeSide` (attachment along an undirected edge).

.. code-block:: python

   # Get all sides (both edge sides and hybrid sides)
   all_sides = generator.sides  # Returns list[Side]
   
   # Get only hybrid sides (hybrid nodes with out-degree 0)
   hybrid_sides = generator.hybrid_sides  # Returns list[HybridSide]
   
   # Get all directed edge sides (all directed edges as DirEdgeSide objects)
   directed_edge_sides = generator.directed_edge_sides  # Returns list[DirEdgeSide]
   
   # Get all undirected edge sides
   undirected_edge_sides = generator.undirected_edge_sides  # Returns list[UndirEdgeSide]
   
   # Get all edge sides (both directed and undirected)
   edge_sides = generator.edge_sides  # Returns list[DirEdgeSide | UndirEdgeSide]
   
   # Get parallel directed edge sides (groups of parallel directed edges)
   parallel_directed_sides = generator.parallel_directed_edge_sides  # Returns list[tuple[DirEdgeSide, ...]]
   
   # Get non-parallel directed edge sides
   non_parallel_directed_sides = generator.non_parallel_directed_edge_sides  # Returns list[DirEdgeSide]

**Validation**

Generators are validated at construction time. The :meth:`validate` method can be called
explicitly to check the generator structure:

.. code-block:: python

   # Validate generator structure
   generator.validate()  # Raises exception if generator is invalid

Generator Functions
-------------------

The generator module provides functions to convert directed generators to semi-directed generators
and generate all possible level-k semi-directed generators.

**Converting Directed Generators**

The :func:`phylozoo.core.network.sdnetwork.generator.dgenerator_to_sdgenerator` function converts
a directed generator to a semi-directed generator by undirecting all non-hybrid edges and
suppressing the degree-2 root node.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import dgenerator_to_sdgenerator
   from phylozoo.core.network.dnetwork.generator import DirectedGenerator
   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   # Create a directed generator
   d_gen_graph = DirectedMultiGraph(edges=[(8, 4), (8, 4)])
   d_gen = DirectedGenerator(d_gen_graph)
   
   # Convert to semi-directed generator
   sd_gen = dgenerator_to_sdgenerator(d_gen)

**Generating All Level-k Generators**

The :func:`phylozoo.core.network.sdnetwork.generator.all_level_k_generators` function generates
all possible level-k semi-directed generators. This is useful for network classification,
as any level-k network can be constructed by combining these generators. The function
returns a set of all distinct generator structures for the specified level.

The function constructs generators by first getting all level-k directed generators and
then converting them to semi-directed generators using `dgenerator_to_sdgenerator`.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import all_level_k_generators
   
   # Get all level-0 generators (single node)
   level_0_generators = all_level_k_generators(0)
   print(f"Number of level-0 generators: {len(level_0_generators)}")
   
   # Get all level-1 generators
   level_1_generators = all_level_k_generators(1)
   print(f"Number of level-1 generators: {len(level_1_generators)}")
   
   # Get all level-2 generators
   level_2_generators = all_level_k_generators(2)
   print(f"Number of level-2 generators: {len(level_2_generators)}")

Generator Sides
---------------

Generator sides represent attachment points where additional structure (leaves, trees) can be
attached when building generators from lower-level generators. The module provides classes
for representing sides, including the base :class:`Side` class and specific side types.

**Side Base Class**

The :class:`phylozoo.core.network.dnetwork.generator.Side` class is the base class for all sides.
It serves as a general wrapper class. Use :class:`HybridSide` for node sides,
:class:`DirEdgeSide` for directed edge sides, and :class:`UndirEdgeSide` for undirected edge sides.

**HybridSide**

The :class:`phylozoo.core.network.dnetwork.generator.HybridSide` class represents a hybrid side of a generator.
A hybrid side is a node with in-degree 2 and out-degree 0, representing a hybrid node that serves
as an attachment point.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import HybridSide
   
   # Create a hybrid side
   hybrid_side = HybridSide(node=3)
   print(hybrid_side.node)  # 3

**DirEdgeSide**

The :class:`phylozoo.core.network.dnetwork.generator.DirEdgeSide` class represents a directed edge side of a generator.
A directed edge side represents a directed edge (possibly parallel) that serves as an attachment point.
The edge is identified by its endpoints and key.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import DirEdgeSide
   
   # Create a directed edge side
   edge_side = DirEdgeSide(u=8, v=4, key=0)
   print(edge_side.u)    # 8
   print(edge_side.v)    # 4
   print(edge_side.key)  # 0

**UndirEdgeSide**

The :class:`phylozoo.core.network.sdnetwork.generator.UndirEdgeSide` class represents an undirected edge side of a generator.
An undirected edge side represents an undirected edge that serves as an attachment point.
The edge is identified by its endpoints and key.

.. code-block:: python

   from phylozoo.core.network.sdnetwork.generator import UndirEdgeSide
   
   # Create an undirected edge side
   undir_edge_side = UndirEdgeSide(u=1, v=2, key=0)
   print(undir_edge_side.u)    # 1
   print(undir_edge_side.v)    # 2
   print(undir_edge_side.key)  # 0

**Working with Sides**

You can iterate over generator sides and check their types:

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator import HybridSide, DirEdgeSide
   from phylozoo.core.network.sdnetwork.generator import UndirEdgeSide
   
   generator = SemiDirectedGenerator(...)
   for side in generator.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, DirEdgeSide):
           print(f"Directed edge side: ({side.u}, {side.v}, key={side.key})")
       elif isinstance(side, UndirEdgeSide):
           print(f"Undirected edge side: ({side.u}, {side.v}, key={side.key})")

.. note::
   Generators are used internally for network classification and construction. They 
   are particularly useful for understanding the structure of level-k networks.

.. tip::
   Use ``all_level_k_generators(k)`` to enumerate all possible level-k generator 
   structures. This is useful for network classification and generation.

See Also
--------

- :doc:`Semi-Directed Network Algorithms <semi_directed_network_algorithms>` - Network level classification
- :doc:`Semi-Directed Network Class <semi_directed_network_class>` - The SemiDirectedPhyNetwork class
- :doc:`Directed Generator <../directed/directed_generator>` - Directed generators
- :doc:`API Reference <../../../../api/core/network/index>` - Complete function signatures and detailed examples
