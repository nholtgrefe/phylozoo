Generators
==========

Generators are minimal biconnected components that represent the core structure of 
level-k phylogenetic networks :cite:`PhyloZoo2024`. They are used to characterize 
and construct networks based on their level.

Directed Generators
-------------------

Directed generators represent the structure of directed level-k networks:

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

Extract generators from networks:

.. code-block:: python

   from phylozoo.core.dnetwork.generator import generators_from_network
   from phylozoo import DirectedPhyNetwork
   
   # Extract generators from a network
   network = DirectedPhyNetwork.load("network.enewick")
   generators = list(generators_from_network(network))
   
   for gen in generators:
       print(f"Level: {gen.level}, Hybrid nodes: {gen.hybrid_nodes}")

Generate all level-k generators:

.. code-block:: python

   # Get all level-2 generators
   level_2_generators = all_level_k_generators(2)
   print(f"Number of level-2 generators: {len(level_2_generators)}")

Semi-Directed Generators
------------------------

Semi-directed generators represent the structure of semi-directed level-k networks:

.. code-block:: python

   from phylozoo.core.sdnetwork.generator import (
       SemiDirectedGenerator, 
       all_level_k_generators,
       dgenerator_to_sdgenerator
   )
   
   # Convert directed generator to semi-directed
   d_gen = DirectedGenerator(...)
   sd_gen = dgenerator_to_sdgenerator(d_gen)
   
   # Get all level-k semi-directed generators
   level_2_sd_generators = all_level_k_generators(2)

Generator Sides
---------------

Generators have "sides" which are attachment points where networks can be connected:

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

Available Classes and Functions
-------------------------------

**Directed Generators:**

* **DirectedGenerator** - Level-k generator for directed phylogenetic networks. 
  Represents minimal biconnected components. See 
  :class:`phylozoo.core.network.dnetwork.generator.DirectedGenerator` for full API.

* **generators_from_network(network)** - Extract generators from a directed network. 
  Network must be binary. Returns iterator of DirectedGenerator objects. See 
  :func:`phylozoo.core.network.dnetwork.generator.generators_from_network`.

* **all_level_k_generators(k)** - Generate all level-k directed generators. Returns 
  set of DirectedGenerator objects. See 
  :func:`phylozoo.core.network.dnetwork.generator.all_level_k_generators`.

**Semi-Directed Generators:**

* **SemiDirectedGenerator** - Level-k generator for semi-directed phylogenetic networks. 
  See :class:`phylozoo.core.network.sdnetwork.generator.SemiDirectedGenerator` for full API.

* **dgenerator_to_sdgenerator(d_generator)** - Convert directed generator to semi-directed. 
  Returns SemiDirectedGenerator. See 
  :func:`phylozoo.core.network.sdnetwork.generator.dgenerator_to_sdgenerator`.

* **all_level_k_generators(k)** - Generate all level-k semi-directed generators. Returns 
  set of SemiDirectedGenerator objects. See 
  :func:`phylozoo.core.network.sdnetwork.generator.all_level_k_generators`.

**Generator Sides:**

* **Side** - Base class for generator sides (attachment points). Abstract base class.

* **HybridSide** - Side representing attachment at a hybrid node. See 
  :class:`phylozoo.core.network.dnetwork.generator.HybridSide`.

* **DirEdgeSide** - Side representing attachment along a directed edge. See 
  :class:`phylozoo.core.network.dnetwork.generator.DirEdgeSide`.

* **UndirEdgeSide** - Side representing attachment along an undirected edge (semi-directed 
  only). See :class:`phylozoo.core.network.sdnetwork.generator.UndirEdgeSide`.

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
   For network level classification, see :doc:`Networks (Advanced) <core_networks_advanced>`. 
   For network construction, see the generator construction functions.
