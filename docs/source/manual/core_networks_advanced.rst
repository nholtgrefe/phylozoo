Networks (Advanced)
===================

This page covers advanced network features, transformations, and analysis capabilities. 
For basic network operations, see :doc:`Networks (Basic) <core_networks_basic>`.

Network Features
----------------

Extract advanced network features:

.. code-block:: python

   from phylozoo.core.dnetwork import features
   
   # Find LSA (Least Stable Ancestor) node
   lsa = features.lsa_node(network)
   
   # Get blobs (maximal biconnected components)
   blob_list = features.blobs(network, trivial=False, leaves=False)
   
   # Get k-blobs (blobs with exactly k incident edges)
   k2_blobs = features.k_blobs(network, k=2)
   
   # Find omnians (nodes that are ancestors of all leaves)
   omnians_set = features.omnians(network)
   
   # Get cut edges and cut vertices
   cut_edges = features.cut_edges(network)
   cut_vertices = features.cut_vertices(network)

For semi-directed networks:

.. code-block:: python

   from phylozoo.core.sdnetwork import features as sd_features
   
   # Get blobs
   blobs = sd_features.blobs(sd_network)
   
   # Get root locations (where network can be rooted)
   root_locs = sd_features.root_locations(sd_network)

Network Classifications
-----------------------

Classify network properties:

.. code-block:: python

   from phylozoo.core.dnetwork import classifications
   
   # Basic properties
   level = classifications.level(network)  # Number of reticulations
   is_binary = classifications.is_binary(network)
   is_tree = classifications.is_tree(network)
   
   # Network types
   is_lsa = classifications.is_lsa_network(network)
   is_galled = classifications.is_galled(network)
   is_treechild = classifications.is_treechild(network)
   is_treebased = classifications.is_treebased(network)
   
   # Check for parallel edges
   has_parallel = classifications.has_parallel_edges(network)
   
   # Ultrametric check
   is_ultrametric = classifications.is_ultrametric(network)

Advanced Transformations
------------------------

Transform networks in various ways:

.. code-block:: python

   from phylozoo.core.dnetwork import transformations
   
   # Convert to LSA network
   lsa_net = transformations.to_lsa_network(network)
   
   # Binary resolution (resolve high-degree nodes)
   binary_net = transformations.binary_resolution(network)
   
   # Identify parallel edges
   identified_net = transformations.identify_parallel_edges(network)
   
   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

Network Derivations
-------------------

Extract derived structures from networks:

.. code-block:: python

   from phylozoo.core.dnetwork import derivations
   
   # Convert to semi-directed network
   sd_net = derivations.to_sd_network(network)
   
   # Extract tree of blobs
   blob_tree = derivations.tree_of_blobs(network)
   
   # Get subnetwork for specific taxa
   subnetwork = derivations.subnetwork(network, taxa={"A", "B", "C"})
   
   # Get k-taxon subnetworks
   k_subnets = derivations.k_taxon_subnetworks(network, k=4)
   
   # Get displayed trees
   for tree in derivations.displayed_trees(network):
       # Process each displayed tree
       pass
   
   # Extract splits and quartets
   splits = derivations.induced_splits(network)
   quartets = derivations.displayed_quartets(network)

Isomorphism Checking
--------------------

Check if two networks have the same structure:

.. code-block:: python

   from phylozoo.core.dnetwork import isomorphism
   
   # Check isomorphism (labels are always checked)
   are_isomorphic = isomorphism.is_isomorphic(net1, net2)
   
   # Check with additional attributes
   are_isomorphic = isomorphism.is_isomorphic(
       net1, net2,
       node_attrs=["custom_attr"],
       edge_attrs=["branch_length"]
   )

Available Functions
-------------------

**Network Features** (``phylozoo.core.dnetwork.features``):

* **lsa_node(network)** - Find the Least Stable Ancestor node. The LSA is the lowest 
  node through which all root-to-leaf paths pass. Returns the node identifier.

* **blobs(network, trivial=True, leaves=True)** - Get all blobs (maximal biconnected 
  components) in the network. Can filter trivial blobs and leaf-only blobs. Returns 
  list of sets of node IDs.

* **k_blobs(network, k)** - Get blobs with exactly k incident edges. Useful for 
  identifying cycles and reticulation structures. Returns list of sets of node IDs.

* **omnians(network)** - Find omnian nodes (nodes that are ancestors of all leaves). 
  Returns set of node IDs.

* **cut_edges(network)** - Get all cut edges (edges whose removal disconnects the 
  network). Returns set of edge tuples.

* **cut_vertices(network)** - Get all cut vertices (nodes whose removal disconnects 
  the network). Returns set of node IDs.

**Network Classifications** (``phylozoo.core.dnetwork.classifications``):

* **level(network)** - Calculate the level (number of reticulations) of the network. 
  Returns integer.

* **is_lsa_network(network)** - Check if network is an LSA network (root equals LSA 
  node). Returns boolean.

* **is_binary(network)** - Check if network is binary (all internal nodes have degree 
  3). Returns boolean.

* **is_tree(network)** - Check if network is actually a tree (no hybrid nodes). 
  Returns boolean.

* **is_galled(network)** - Check if network is galled (each reticulation is in its own 
  cycle). Returns boolean.

* **is_treechild(network)** - Check if network is tree-child (each internal node has 
  at least one tree child). Returns boolean.

* **is_treebased(network)** - Check if network is tree-based (can be obtained from a 
  tree by adding edges). Returns boolean.

* **is_ultrametric(network)** - Check if network is ultrametric (all root-to-leaf 
  distances equal). Returns boolean.

* **has_parallel_edges(network)** - Check if network contains parallel edges. Returns 
  boolean.

**Network Transformations** (``phylozoo.core.dnetwork.transformations``):

* **to_lsa_network(network)** - Convert network to LSA form by suppressing nodes above 
  the LSA. Returns new DirectedPhyNetwork.

* **binary_resolution(network)** - Resolve high-degree nodes to binary form using 
  caterpillar structures. Preserves gamma values and branch lengths. Returns new 
  DirectedPhyNetwork.

* **identify_parallel_edges(network)** - Identify parallel edges by adding explicit 
  edge keys. Returns new DirectedPhyNetwork.

* **suppress_2_blobs(network)** - Suppress degree-2 nodes in 2-blobs. Simplifies 
  network structure. Returns new DirectedPhyNetwork.

**Network Derivations** (``phylozoo.core.dnetwork.derivations``):

* **to_sd_network(network)** - Convert directed network to semi-directed network. 
  Undirects non-hybrid edges. Returns SemiDirectedPhyNetwork.

* **tree_of_blobs(network)** - Extract tree structure of blobs. Returns 
  DirectedPhyNetwork representing blob tree.

* **subnetwork(network, taxa)** - Extract subnetwork induced by set of taxa. Returns 
  new DirectedPhyNetwork.

* **k_taxon_subnetworks(network, k)** - Get all k-taxon subnetworks. Returns iterator 
  of DirectedPhyNetwork objects.

* **displayed_trees(network, probability=False)** - Get all displayed trees (trees 
  embedded in network). If probability=True, includes edge probabilities. Returns 
  iterator of DirectedPhyNetwork objects.

* **induced_splits(network)** - Extract splits induced by network. Returns SplitSystem.

* **displayed_quartets(network)** - Extract quartet profiles from displayed trees. 
  Returns QuartetProfileSet.

**Isomorphism** (``phylozoo.core.dnetwork.isomorphism``):

* **is_isomorphic(net1, net2, node_attrs=None, edge_attrs=None, graph_attrs=None)** - 
  Check if two networks are isomorphic. Labels are always checked. Additional 
  attributes can be specified. Returns boolean.

**Semi-Directed Network Features** (``phylozoo.core.sdnetwork.features``):

* **blobs(network, trivial=True, leaves=True)** - Get blobs in semi-directed network. 
  Returns list of sets of node IDs.

* **root_locations(network)** - Find all possible root locations (nodes or edges). 
  Returns list of root locations.

.. note::
   All transformation and derivation functions return new network instances. Networks 
   are immutable - see :doc:`Advanced Topics <advanced>` for details.

.. tip::
   Use network classifications to determine which algorithms can be applied to your 
   network. For example, certain diversity optimization methods require networks 
   without 2-blobs or parallel edges.

.. warning::
   Some advanced operations (like binary resolution) can significantly increase network 
   size. Use with caution on large networks.
