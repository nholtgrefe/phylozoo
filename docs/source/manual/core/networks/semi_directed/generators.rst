Semi-Directed Network Generators
==================================

The :mod:`phylozoo.core.network.sdnetwork.generator` module provides the :class:`SemiDirectedGenerator`
class and related functions for working with level-k semi-directed network generators. Generators
are minimal biconnected components that represent the core structure of level-k semi-directed
phylogenetic networks :cite:`PhyloZoo2024`. They are used to characterize and construct networks
based on their level.

All classes and functions on this page can be imported from the core network generator module:

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import SemiDirectedGenerator, all_level_k_generators, generators_from_network, dgenerator_to_sdgenerator
   # or directly
   from phylozoo.core.sdnetwork.generator import SemiDirectedGenerator

Working with Semi-Directed Generators
---------------------------------------

Semi-directed generators represent the fundamental building blocks of level-k semi-directed
phylogenetic networks. Each generator is a minimal biconnected component that can be combined
with other generators to construct networks of a specific level.

Creating Generators
^^^^^^^^^^^^^^^^^^^

Semi-directed generators are minimal biconnected components that represent the core
structure of level-k semi-directed phylogenetic networks. They can be created by
converting directed generators using :func:`dgenerator_to_sdgenerator`, or extracted
directly from semi-directed networks. Each generator has a level (k), hybrid nodes,
and sides (attachment points) where networks can be connected.

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import (
       SemiDirectedGenerator, 
       all_level_k_generators,
       dgenerator_to_sdgenerator
   )
   from phylozoo.core.dnetwork.generator import DirectedGenerator
   
   # Convert directed generator to semi-directed
   d_gen = DirectedGenerator(...)
   sd_gen = dgenerator_to_sdgenerator(d_gen)
   
   # Access properties
   level = sd_gen.level
   hybrid_nodes = sd_gen.hybrid_nodes
   sides = sd_gen.sides  # List of attachment points

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.sdnetwork.generator.generators_from_network` function extracts
all generators from a binary semi-directed network. This function decomposes the network
into its biconnected components and identifies the generator structures. The network
must be binary (no parallel edges or high-degree nodes) for accurate generator extraction.

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import generators_from_network
   from phylozoo import SemiDirectedPhyNetwork
   
   # Extract generators from a network
   network = SemiDirectedPhyNetwork.load("network.newick")
   generators = list(generators_from_network(network))
   
   for gen in generators:
       print(f"Level: {gen.level}, Hybrid nodes: {gen.hybrid_nodes}")

Generating All Level-k Generators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.sdnetwork.generator.all_level_k_generators` function generates
all possible level-k semi-directed generators. This is useful for network classification,
as any level-k network can be constructed by combining these generators. The function
returns a set of all distinct generator structures for the specified level.

.. code-block:: python

   # Get all level-2 semi-directed generators
   level_2_generators = all_level_k_generators(2)
   print(f"Number of level-2 generators: {len(level_2_generators)}")

Generator Sides
^^^^^^^^^^^^^^^

Generators have "sides" which are attachment points where networks can be connected
to form larger networks. There are two types of sides: HybridSide (attachment at a
hybrid node) and UndirEdgeSide (attachment along an undirected edge). Sides define
how generators can be combined to construct level-k networks.

**HybridSide**

The :class:`phylozoo.core.sdnetwork.generator.HybridSide` class represents attachment
at a hybrid node.

**UndirEdgeSide**

The :class:`phylozoo.core.sdnetwork.generator.UndirEdgeSide` class represents attachment
along an undirected edge.

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import Side, HybridSide, UndirEdgeSide
   
   generator = SemiDirectedGenerator(...)
   for side in generator.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, UndirEdgeSide):
           print(f"Undirected edge side: {side.edge}")

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

- :doc:`Semi-Directed Networks (Advanced) <advanced>` - Network level classification
- :doc:`Directed Network Generators <../directed/generators>` - Directed generators
