Semi-Directed Network Generators
==================================

Generators are minimal biconnected components that represent the core structure of 
level-k semi-directed phylogenetic networks :cite:`PhyloZoo2024`. They are used to characterize 
and construct networks based on their level.

Creating Generators
-------------------

Semi-directed generators represent the structure of semi-directed level-k networks:

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

Extracting Generators from Networks
------------------------------------

Extract generators from networks:

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import generators_from_network
   from phylozoo import SemiDirectedPhyNetwork
   
   # Extract generators from a network
   network = SemiDirectedPhyNetwork.load("network.newick")
   generators = list(generators_from_network(network))
   
   for gen in generators:
       print(f"Level: {gen.level}, Hybrid nodes: {gen.hybrid_nodes}")

Generating All Level-k Generators
---------------------------------

Generate all level-k generators:

.. code-block:: python

   # Get all level-2 semi-directed generators
   level_2_generators = all_level_k_generators(2)
   print(f"Number of level-2 generators: {len(level_2_generators)}")

Generator Sides
---------------

Generators have "sides" which are attachment points where networks can be connected:

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import Side, HybridSide, UndirEdgeSide
   
   # Sides represent attachment points
   # HybridSide: attachment at a hybrid node
   # UndirEdgeSide: attachment along an undirected edge
   
   generator = SemiDirectedGenerator(...)
   for side in generator.sides:
       if isinstance(side, HybridSide):
           print(f"Hybrid side at node: {side.node}")
       elif isinstance(side, UndirEdgeSide):
           print(f"Undirected edge side: {side.edge}")

API Reference
-------------

**Classes:**

* **SemiDirectedGenerator** - Level-k generator for semi-directed phylogenetic networks. 
  Represents minimal biconnected components. See 
  :class:`phylozoo.core.network.sdnetwork.generator.SemiDirectedGenerator` for full API.

**Functions:**

* **dgenerator_to_sdgenerator(d_generator)** - Convert directed generator to semi-directed. 
  Returns SemiDirectedGenerator. See 
  :func:`phylozoo.core.network.sdnetwork.generator.dgenerator_to_sdgenerator`.

* **generators_from_network(network)** - Extract generators from a semi-directed network. 
  Network must be binary. Returns iterator of SemiDirectedGenerator objects. See 
  :func:`phylozoo.core.network.sdnetwork.generator.generators_from_network`.

* **all_level_k_generators(k)** - Generate all level-k semi-directed generators. Returns 
  set of SemiDirectedGenerator objects. See 
  :func:`phylozoo.core.network.sdnetwork.generator.all_level_k_generators`.

**Generator Sides:**

* **Side** - Base class for generator sides (attachment points). Abstract base class.

* **HybridSide** - Side representing attachment at a hybrid node. See 
  :class:`phylozoo.core.network.sdnetwork.generator.HybridSide`.

* **UndirEdgeSide** - Side representing attachment along an undirected edge. See 
  :class:`phylozoo.core.network.sdnetwork.generator.UndirEdgeSide`.

**Generator Properties:**

* **level** - The level (k) of the generator. Returns integer.

* **hybrid_nodes** - Set of hybrid nodes in the generator. Returns set of node IDs.

* **sides** - List of sides (attachment points). Returns list of Side objects.

.. note::
   Generators are used internally for network classification and construction. They 
   are particularly useful for understanding the structure of level-k networks.

.. tip::
   Use ``all_level_k_generators(k)`` to enumerate all possible level-k generator 
   structures. This is useful for network classification and generation.

.. warning::
   Extracting generators from networks requires the network to be binary. Networks 
   with parallel edges may produce unexpected results.

.. seealso::
   For network level classification, see :doc:`Semi-Directed Networks (Advanced) <advanced>`. 
   For directed generators, see :doc:`Directed Network Generators <../directed/generators>`.
