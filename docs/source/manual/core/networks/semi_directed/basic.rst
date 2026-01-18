Semi-Directed Networks (Basic)
================================

This page covers basic operations for working with **SemiDirectedPhyNetwork**. For advanced 
features like network analysis, transformations, and classifications, see 
:doc:`Semi-Directed Networks (Advanced) <advanced>`.

SemiDirectedPhyNetwork
-----------------------

**SemiDirectedPhyNetwork** represents networks with directed hybrid edges and undirected 
tree edges. These allow for more flexible representation and are useful for unrooted 
analyses. Semi-directed networks cannot have undirected cycles.

**I/O Formats**: Newick (default, extensions: ``.nwk``, ``.newick``, ``.enewick``, ``.eNewick``, ``.enw``), 
PhyloZoo-DOT (extension: ``.pzdot``). See :doc:`I/O <../../../../io>` for details.

.. figure:: ../../../../_static/images/example_tree_semidirected.png
   :alt: Example semi-directed tree
   :width: 400px
   :align: center
   
   Example of a simple semi-directed tree.

.. figure:: ../../../../_static/images/example_hybrid_semidirected.png
   :alt: Example semi-directed network with hybrid
   :width: 400px
   :align: center
   
   Example of a semi-directed network with a hybrid node.

Creating Networks
-----------------

Create networks from directed and undirected edges:

.. code-block:: python

   from phylozoo import SemiDirectedPhyNetwork
   
   # Simple tree (all edges undirected)
   network = SemiDirectedPhyNetwork(
       undirected_edges=[(3, 1), (3, 2), (3, 4)],
       nodes=[
           (1, {"label": "A"}),
           (2, {"label": "B"}),
           (4, {"label": "C"})
       ]
   )

Create networks with hybrid nodes:

.. code-block:: python

   # Network with hybridization
   hybrid_net = SemiDirectedPhyNetwork(
       directed_edges=[
           (5, 4, {"gamma": 0.6}),  # Hybrid edge
           (6, 4, {"gamma": 0.4})   # Hybrid edge (must sum to 1.0)
       ],
       undirected_edges=[(4, 1), (4, 2), (4, 3)],  # Tree edges
       nodes=[
           (1, {"label": "A"}),
           (2, {"label": "B"}),
           (3, {"label": "C"})
       ]
   )

Loading and Saving
------------------

Networks can be loaded from and saved to files:

.. code-block:: python

   # Save network
   network.save("network.newick")
   
   # Load network
   network = SemiDirectedPhyNetwork.load("network.newick")

See the :doc:`I/O <../../../../io>` page for supported formats and detailed I/O operations.

Basic Properties
----------------

Access basic network properties:

.. code-block:: python

   # Node and edge counts
   num_nodes = network.num_nodes
   num_edges = network.num_edges
   
   # Network structure
   leaves = network.leaves  # Set of leaf node IDs
   taxa = network.taxa       # Set of taxon labels (strings)
   internal_nodes = network.internal_nodes
   
   # Node types
   hybrid_nodes = network.hybrid_nodes
   tree_nodes = network.tree_nodes
   
   # Check if tree
   is_tree = network.is_tree()
   
   # Get root locations (where network can be rooted)
   from phylozoo.core.sdnetwork import features
   root_locs = features.root_locations(network)

Accessing Edge Attributes
-------------------------

Networks support edge attributes like branch lengths, bootstrap values, and gamma 
probabilities (for hybrid edges only):

.. code-block:: python

   # Get branch length (works for both directed and undirected edges)
   bl = network.get_branch_length(3, 1)
   
   # Get bootstrap support
   bs = network.get_bootstrap(3, 1)
   
   # Get gamma (for hybrid edges only)
   gamma = network.get_gamma(5, 4)

Simple Transformations
----------------------

Basic network transformations:

.. code-block:: python

   from phylozoo.core.sdnetwork.transformations import suppress_2_blobs
   
   # Suppress degree-2 nodes (simplify network)
   simplified = suppress_2_blobs(network)

For more advanced transformations, see :doc:`Semi-Directed Networks (Advanced) <advanced>`.

API Reference
-------------

**Class**: :class:`phylozoo.core.network.sdnetwork.SemiDirectedPhyNetwork`

**Basic Properties:**

* **num_nodes** - Number of nodes in the network
* **num_edges** - Number of edges in the network
* **leaves** - Set of leaf node IDs
* **taxa** - Set of taxon labels (strings)
* **internal_nodes** - Set of internal node IDs
* **hybrid_nodes** - Set of hybrid node IDs
* **tree_nodes** - Set of tree node IDs

**Basic Methods:**

* **is_tree()** - Check if network is a tree
* **get_branch_length(u, v, key=None)** - Get branch length for edge
* **get_bootstrap(u, v, key=None)** - Get bootstrap value for edge
* **get_gamma(u, v, key=None)** - Get gamma probability for hybrid edge
* **save(filename)** - Save network to file
* **load(filename)** - Load network from file (class method)

**Basic Transformations:**

* **suppress_2_blobs(network)** - Suppress degree-2 nodes in 2-blobs. Simplifies network 
  structure by removing unnecessary degree-2 nodes while preserving topology. See 
  :func:`phylozoo.core.network.sdnetwork.transformations.suppress_2_blobs` for details.

.. note::
   For advanced network features like blobs, network classifications, and root locations, 
   see :doc:`Semi-Directed Networks (Advanced) <advanced>`.

.. tip::
   All networks are immutable. To modify a network, create a new instance with the 
   desired changes.

.. seealso::
   For directed networks, see :doc:`Directed Networks <../directed/overview>`. 
   For I/O operations, see :doc:`I/O <../../../../io>`.
