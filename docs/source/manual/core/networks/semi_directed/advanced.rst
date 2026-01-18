Semi-Directed Networks (Advanced)
==================================

This page covers advanced features, transformations, and analysis capabilities for 
**SemiDirectedPhyNetwork**. For basic network operations, see :doc:`Semi-Directed Networks (Basic) <basic>`.

Network Features
----------------

Extract advanced network features:

.. code-block:: python

   from phylozoo.core.sdnetwork import features
   
   # Get blobs (maximal biconnected components)
   blob_list = features.blobs(network, trivial=False, leaves=False)
   
   # Get root locations (where network can be rooted)
   root_locs = features.root_locations(network)
   
   # Get up-down paths
   from phylozoo.core.sdnetwork.m_multigraph.features import updown_path_vertices
   path = updown_path_vertices(network._graph, u, v)

Network Classifications
------------------------

Classify network properties:

.. code-block:: python

   from phylozoo.core.sdnetwork import classifications
   
   # Basic properties
   level = classifications.level(network)  # Number of reticulations
   is_binary = classifications.is_binary(network)
   is_tree = classifications.is_tree(network)
   
   # Network types
   is_galled = classifications.is_galled(network)

Advanced Transformations
------------------------

Transform networks in various ways:

.. code-block:: python

   from phylozoo.core.sdnetwork import transformations
   
   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

Network Derivations
-------------------

Extract derived structures from networks:

.. code-block:: python

   from phylozoo.core.sdnetwork import derivations
   
   # Extract tree of blobs
   blob_tree = derivations.tree_of_blobs(network)
   
   # Get subnetwork for specific taxa
   subnetwork = derivations.subnetwork(network, taxa={"A", "B", "C"})
   
   # Get k-taxon subnetworks
   k_subnets = derivations.k_taxon_subnetworks(network, k=4)

API Reference
-------------

**Network Features** (``phylozoo.core.sdnetwork.features``):

* **blobs(network, trivial=True, leaves=True)** - Get all blobs (maximal biconnected 
  components) in the network. Can filter trivial blobs and leaf-only blobs. Returns 
  list of sets of node IDs.

* **root_locations(network)** - Find all possible root locations (nodes or edges). 
  Returns list of root locations.

**Network Classifications** (``phylozoo.core.sdnetwork.classifications``):

* **level(network)** - Calculate the level (number of reticulations) of the network. 
  Returns integer.

* **is_binary(network)** - Check if network is binary (all internal nodes have degree 
  3). Returns boolean.

* **is_tree(network)** - Check if network is actually a tree (no hybrid nodes). 
  Returns boolean.

* **is_galled(network)** - Check if network is galled (no hybrid node is ancestral to 
  another hybrid node in the same blob). Returns boolean.

**Network Transformations** (``phylozoo.core.sdnetwork.transformations``):

* **suppress_2_blobs(network)** - Suppress degree-2 nodes in 2-blobs. Simplifies 
  network structure. Returns new SemiDirectedPhyNetwork.

**Network Derivations** (``phylozoo.core.sdnetwork.derivations``):

* **tree_of_blobs(network)** - Extract tree structure of blobs. Returns 
  SemiDirectedPhyNetwork representing blob tree.

* **subnetwork(network, taxa)** - Extract subnetwork induced by set of taxa. Returns 
  new SemiDirectedPhyNetwork.

* **k_taxon_subnetworks(network, k)** - Get all k-taxon subnetworks. Returns iterator 
  of SemiDirectedPhyNetwork objects.

.. note::
   All transformation and derivation functions return new network instances. Networks 
   are immutable.

.. tip::
   Use network classifications to determine which algorithms can be applied to your 
   network.

.. seealso::
   For basic operations, see :doc:`Semi-Directed Networks (Basic) <basic>`. 
   For directed networks, see :doc:`Directed Networks <../directed/overview>`.
