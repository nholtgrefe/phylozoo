Directed Network (Advanced Features)
====================================

The :mod:`phylozoo.core.network.dnetwork` module provides advanced features, transformations,
and analysis capabilities for :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`. This page covers network features,
classifications, transformations, derivations, and isomorphism checking. For the DirectedPhyNetwork
class and its properties, see :doc:`Directed Network Class <directed_network_class>`.

Structural Features
-------------------

The features module provides functions to extract advanced structural properties of
directed networks. These include LSA (Least Stable Ancestor) nodes, blobs (maximal
biconnected components), omnian nodes, and cut edges/vertices that characterize
network topology.


**LSA Node**

The :func:`~phylozoo.core.network.dnetwork.features.lsa_node` function finds the Least Stable
Ancestor node, which is the lowest node through which all root-to-leaf paths pass.
The LSA is a fundamental concept in phylogenetic network analysis, representing the
most recent common ancestor of all leaves. For convenience, this function is also available as 
the cached :attr:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork.LSA_node` property on the DirectedPhyNetwork class.

.. code-block:: python

   from phylozoo.core.network.dnetwork import features
   
   # Find LSA (Least Stable Ancestor) node
   lsa = features.lsa_node(network)

**Blobs**

The :func:`~phylozoo.core.network.dnetwork.features.blobs` function extracts all blobs (maximal
biconnected components) in the network, represented by sets of nodes. Blobs represent the reticulate structure of
the network, where cycles and hybrid nodes create biconnected regions. The function can optionally exclude trivial (single-node) blobs and/or blobs that contain only leaves.

.. code-block:: python

   # Get all blobs (including trivial and leaf blobs)
   all_blobs = features.blobs(network)
   
   # Get only non-trivial blobs (excluding single-node blobs)
   non_trivial_blobs = features.blobs(network, trivial=False)
   
   # Get blobs excluding leaf-only blobs
   internal_blobs = features.blobs(network, leaves=False)

**k-Blobs**

The :func:`~phylozoo.core.network.dnetwork.features.k_blobs` function gets blobs with exactly
:math:`k` incident edges, represented by sets of nodes. An incident edge is any edge that
connects a node inside the blob to a node outside the blob.

.. code-block:: python

   # Get 2-blobs (blobs with exactly 2 incident edges)
   k2_blobs = features.k_blobs(network, k=2)

**Omnians**

The :func:`~phylozoo.core.network.dnetwork.features.omnians` function finds omnian nodes,
which are nodes with every child being a hybrid node :cite:`Jetten2016`.

.. code-block:: python

   # Find omnians (nodes that are ancestors of all leaves)
   omnians_set = features.omnians(network)

**Cut Edges and Vertices**

The :func:`~phylozoo.core.network.dnetwork.features.cut_edges` and :func:`~phylozoo.core.network.dnetwork.features.cut_vertices`
functions identify edges and vertices whose removal disconnects the network. These
are important for understanding network connectivity and decomposition.

.. code-block:: python

   # Get cut edges and cut vertices
   cut_edges = features.cut_edges(network)
   cut_vertices = features.cut_vertices(network)

Classifications
-----------------------

The classifications module provides functions to determine various mathematical and
phylogenetic properties of networks. These include level (number of reticulations),
network types (LSA, galled, tree-child, tree-based, normal), and structural properties
like binary resolution and ultrametricity.

Basic Network Properties
^^^^^^^^^^^^^^^^^^^^^^^^^

**Phylogenetic Trees**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_tree` function checks if the network
is a rooted phylogenetic tree (no hybrid nodes). For convenience, this function is also available as a class method on the DirectedPhyNetwork class:

.. code-block:: python

   # Check if network is a tree
   is_tree = network.is_tree()
   # or
   is_tree = classifications.is_tree(network)

**Binary Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_binary` function checks if the network
is binary (all internal nodes have degree 3, except for the root node which must have degree 2).

.. code-block:: python

   # Check if network is binary
   is_binary = classifications.is_binary(network)

**Simple Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_simple` function checks if the network
is simple: containing exactly one non-leaf blob.

.. code-block:: python

   # Check if network is simple
   is_simple = classifications.is_simple(network)

**LSA Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_lsa_network` function checks if the
network is an LSA network, meaning the root equals the LSA node. This is a desirable
property for many network analyses.

.. code-block:: python

   # Check if network is LSA network
   is_lsa = classifications.is_lsa_network(network)

**Parallel Edges**

The :func:`~phylozoo.core.network.dnetwork.classifications.has_parallel_edges` function checks if
the network contains parallel edges (multiple edges between the same pair of nodes).

.. code-block:: python

   # Check for parallel edges
   has_parallel = classifications.has_parallel_edges(network)


Network Classes
^^^^^^^^^^^^^^^

**Reticulation Number**

The :func:`~phylozoo.core.network.dnetwork.classifications.reticulation_number` function calculates
the reticulation number, which is the total number of hybrid edges minus the total number
of hybrid nodes.

.. code-block:: python

   # Calculate reticulation number
   reticulation_number = classifications.reticulation_number(network)

**Level and Vertex Level**

The :func:`~phylozoo.core.network.dnetwork.classifications.level` and :func:`~phylozoo.core.network.dnetwork.classifications.vertex_level` functions calculate the level and vertex level of the network, respectively.

The level is the maximum over all blobs of (number of hybrid edges minus number of hybrid nodes) in that blob.
The vertex level is the maximum over all blobs of the number of hybrid nodes in that blob.
Note that the vertex level is always less than or equal to the level and they coincide for binary networks.

.. code-block:: python

   # Calculate level and vertex level
   level = classifications.level(network)
   vertex_level = classifications.vertex_level(network)

**Stack-Free Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_stackfree` function checks if the network
is stack-free, meaning no two hybrid nodes share a common parent.

.. code-block:: python

   # Check if network is stack-free
   is_stackfree = classifications.is_stackfree(network)

**Tree-Child Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_treechild` function checks if the
network is tree-child, meaning each internal node has at least one tree child. This
property ensures that the network has a tree-like structure.

.. code-block:: python

   # Check if network is tree-child
   is_treechild = classifications.is_treechild(network)

**Tree-Based and Strictly Tree-based Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_treebased` and :func:`~phylozoo.core.network.dnetwork.classifications.is_strictly_treebased` functions check if the
network is tree-based and strictly tree-based. This is currently a stub function and will be implemented in the future.

.. code-block:: python

   # Check if network is tree-based
   is_treebased = classifications.is_treebased(network)
   is_strictly_treebased = classifications.is_strictly_treebased(network)

**Galled Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_galled` function checks if the network
is galled, meaning each reticulation is in its own cycle. Galled networks have special
properties that make them easier to analyze.

.. code-block:: python

   # Check if network is galled
   is_galled = classifications.is_galled(network)

**Normal Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_normal` function checks if the network
is normal, meaning it is tree-child without shortcuts. Normal networks have additional
structural constraints.

.. code-block:: python

   # Check if network is normal
   is_normal = classifications.is_normal(network)

**Ultrametric Networks**

The :func:`~phylozoo.core.network.dnetwork.classifications.is_ultrametric` function checks if all
root-to-leaf distances are equal, using the branch lengths of the edges.

.. code-block:: python

   # Check if network is ultrametric
   is_ultrametric = classifications.is_ultrametric(network)

Transformations
------------------------

The transformations module provides functions to modify network structures. These include converting to LSA form,
binary resolution of high-degree nodes, identifying parallel edges, and suppressing 2-blobs.

**LSA Network Conversion**

The :func:`~phylozoo.core.network.dnetwork.transformations.to_lsa_network` function converts
a network to LSA form by removing all nodes above the LSA.

.. code-block:: python

   from phylozoo.core.network.dnetwork import transformations
   
   # Convert to LSA network
   lsa_net = transformations.to_lsa_network(network)

**Parallel Edge Identification**

The :func:`~phylozoo.core.network.dnetwork.transformations.identify_parallel_edges` function
identifies parallel edges by replacing all parallel edges between the same pair of nodes with a single edge, and suppressing all resulting degree-2 nodes.

.. code-block:: python

   # Identify parallel edges
   identified_net = transformations.identify_parallel_edges(network)

**Suppressing 2-Blobs**

The :func:`~phylozoo.core.network.dnetwork.transformations.suppress_2_blobs` function suppresses
all 2-blobs and the resulting degree-2 nodes.

.. code-block:: python

   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

**Binary Resolution**

The :func:`~phylozoo.core.network.dnetwork.transformations.binary_resolution` function resolves
high-degree nodes to binary form using caterpillar structures, preserving gamma values
and branch lengths. This transformation is useful for algorithms that require binary networks.

.. code-block:: python

   # Binary resolution (resolve high-degree nodes)
   binary_net = transformations.binary_resolution(network)


Derivations
-------------------

The derivations module provides functions to extract derived structures from networks,
including conversion to semi-directed networks, tree-of-blobs, subnetworks for specific
taxa sets, displayed trees, and induced splits and quartets.

Phylogenetic Trees
^^^^^^^^^^^^^^^^^^^

**Tree-of-Blobs**

The :func:`~phylozoo.core.network.dnetwork.derivations.tree_of_blobs` function extracts the
tree structure of blobs, representing the high-level topology of the network. This
simplifies the network by collapsing each blob into a single node.

.. code-block:: python

   # Extract tree of blobs
   tob = derivations.tree_of_blobs(network)

**Displayed Trees**

The :func:`~phylozoo.core.network.dnetwork.derivations.displayed_trees` function generates all
displayed trees (trees embedded in the network). A displayed tree is obtained by taking
a switching (deleting all but one parent edge per hybrid node), then removing degree-1 nodes
and suppressing degree-2 nodes. Optionally the probability of each displayed tree is stored in the network's 'probability' attribute.

.. code-block:: python

   # Get displayed trees without probabilities
   for tree in derivations.displayed_trees(network):
       # Process each displayed tree
       pass
   
   # Get displayed trees with probabilities
   for tree in derivations.displayed_trees(network, probability=True):
       prob = tree.get_network_attribute('probability')
       print(f"Tree probability: {prob}")

Phylogenetic Networks
^^^^^^^^^^^^^^^^^^^^^^

**Subnetworks**

The :func:`~phylozoo.core.network.dnetwork.derivations.subnetwork` function extracts a subnetwork
induced by a specific set of taxa. The subnetwork is defined as the union of all directed
paths from the requested leaves up to the root (i.e., all their ancestors and the leaves
themselves). The induced subgraph is taken, then degree-2 internal nodes are suppressed.

The function has three post-processing options: it can suppress 2-blobs, identify/merge parallel edges, and/or convert the result to an LSA-network (with the LSA as root).

.. code-block:: python

   # Get basic subnetwork
   subnetwork = derivations.subnetwork(network, taxa=["A", "B", "C"])
   
**k-Taxon Subnetworks**

The :func:`~phylozoo.core.network.dnetwork.derivations.k_taxon_subnetworks` function generates all
subnetworks induced by exactly k taxa. For each combination of k taxa, the corresponding
subnetwork is computed using the `subnetwork` function.

.. code-block:: python

   # Get all 4-taxon subnetworks
   k_subnets = derivations.k_taxon_subnetworks(network, k=4)

**Conversion to Semi-Directed Network**

The :func:`~phylozoo.core.network.dnetwork.derivations.to_sd_network` function converts a
directed network to a :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` by undirecting non-hybrid edges and suppressing the root node.

.. code-block:: python

   from phylozoo.core.network.dnetwork import derivations
   
   # Convert to semi-directed network
   sd_net = derivations.to_sd_network(network)

Pairwise Distances
^^^^^^^^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.network.dnetwork.derivations.distances` function computes pairwise distances
between taxa. The function can compute shortest, longest, or
probability-weighted average distances across all displayed trees of the network. This returns a 
:class:`~phylozoo.core.distance.base.DistanceMatrix` object.

.. code-block:: python

   # Compute distance matrix with average distances
   from phylozoo.core.network.dnetwork.derivations import distances
   distance_matrix = distances(network, mode='average')  # or 'shortest', 'longest'

Partitions and Splits
^^^^^^^^^^^^^^^^^^^^^

**Partition from Blob**

The :func:`~phylozoo.core.network.dnetwork.derivations.partition_from_blob` function extracts the partition
of taxa induced by removing a specific blob from the network. 

.. code-block:: python

   # Extract partition from specific blob
   from phylozoo.core.network.dnetwork.features import blobs
   blob = list(blobs(network))[0]
   partition = derivations.partition_from_blob(network, blob)

**Split from Cut-Edge**

The :func:`~phylozoo.core.network.dnetwork.derivations.split_from_cutedge` function extracts the split
induced by a specific cut edge. The split is the 2-partition of taxa obtained when removing the edge from the network.

.. code-block:: python

   # Extract split from specific cut edge
   from phylozoo.core.network.dnetwork.features import cut_edges
   cut_edge = list(cut_edges(network))[0]
   split = derivations.split_from_cutedge(network, cut_edge)

**Induced Splits**

The :func:`~phylozoo.core.network.dnetwork.derivations.induced_splits` function extracts all splits
induced by cut-edges of the network. This returns a :class:`~phylozoo.core.split.splitsystem.SplitSystem` object.

.. code-block:: python

   # Extract induced splits
   splits = derivations.induced_splits(network)

**Displayed Splits**

The :func:`~phylozoo.core.network.dnetwork.derivations.displayed_splits` function extracts all splits
induced by all displayed trees, weighted by their probabilities. This returns a :class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem` object.

.. code-block:: python

   # Extract displayed splits
   displayed_splits = derivations.displayed_splits(network)

Quartets
^^^^^^^^

The :func:`~phylozoo.core.network.dnetwork.derivations.displayed_quartets` function extracts quartet
profiles from displayed trees. This returns a :class:`~phylozoo.core.quartet.qprofileset.QuartetProfileSet` object.

.. code-block:: python

   # Extract displayed quartets
   quartets = derivations.displayed_quartets(network)


Isomorphism Checking
--------------------

The isomorphism module provides functions to check if two networks have the same
topological structure. 

The :func:`~phylozoo.core.network.dnetwork.isomorphism.is_isomorphic` function checks if two
networks are isomorphic. Labels are always checked, and additional attributes can
be specified for comparison.

.. code-block:: python

   from phylozoo.core.network.dnetwork import isomorphism
   
   # Check isomorphism (labels are always checked)
   are_isomorphic = isomorphism.is_isomorphic(net1, net2)
   
   # Check with additional attributes
   are_isomorphic = isomorphism.is_isomorphic(
       net1, net2,
       node_attrs=["custom_attr"],
       edge_attrs=["branch_length"]
   )


See Also
--------

- :doc:`API Reference <../../../../api/core/network/dnetwork/index>` - Complete function signatures and detailed examples
- :doc:`Directed Network Class <directed_network_class>` - The DirectedPhyNetwork class and its properties
- :doc:`Directed Generator <directed_generator>` - Level-k network generators
- :doc:`Semi-Directed Networks <../semi_directed/overview>` - Semi-directed network representations
