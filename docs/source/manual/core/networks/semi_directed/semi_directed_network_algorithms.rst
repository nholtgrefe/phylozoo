Semi-Directed Network (Advanced Features)
==========================================

The :mod:`phylozoo.core.network.sdnetwork` module provides advanced features, transformations,
and analysis capabilities for :class:`SemiDirectedPhyNetwork`. This page covers network features,
classifications, transformations, derivations, and isomorphism checking. For the SemiDirectedPhyNetwork
class and its properties, see :doc:`Semi-Directed Network Class <semi_directed_network_class>`.

Structural Features
-------------------

The features module provides functions to extract advanced structural properties of
semi-directed networks. These include blobs (maximal biconnected components), root
locations (where the network can be rooted), omnian nodes, and cut edges/vertices that characterize
network topology.

**Blobs**

The :func:`phylozoo.core.network.sdnetwork.features.blobs` function extracts all blobs (maximal
biconnected components) in the network, represented by sets of nodes. Blobs represent the reticulate structure of
the network, where cycles and hybrid nodes create biconnected regions. The function can optionally exclude trivial (single-node) blobs and/or blobs that contain only leaves.

.. code-block:: python

   from phylozoo.core.network.sdnetwork import features
   
   # Get all blobs (including trivial and leaf blobs)
   all_blobs = features.blobs(network)
   
   # Get only non-trivial blobs (excluding single-node blobs)
   non_trivial_blobs = features.blobs(network, trivial=False)
   
   # Get blobs excluding leaf-only blobs
   internal_blobs = features.blobs(network, leaves=False)

**k-Blobs**

The :func:`phylozoo.core.network.sdnetwork.features.k_blobs` function gets blobs with exactly
:math:`k` incident edges, represented by sets of nodes. An incident edge is any edge that
connects a node inside the blob to a node outside the blob.

.. code-block:: python

   # Get k-blobs (blobs with exactly k incident edges)
   k2_blobs = features.k_blobs(network, k=2)

**Root Locations**

The :func:`phylozoo.core.network.sdnetwork.features.root_locations` function finds all possible
root locations (nodes or edges) where the network can be rooted. Unlike directed networks,
semi-directed networks do not have a fixed root, but can be rooted at various locations.

.. code-block:: python

   # Get root locations (where network can be rooted)
   root_locs = features.root_locations(network)

**Omnians**

The :func:`phylozoo.core.network.sdnetwork.features.omnians` function finds omnian nodes,
which are nodes with every child being a hybrid node :cite:`jetten2016nonbinary`.

.. code-block:: python

   # Find omnians (nodes that are ancestors of all leaves)
   omnians_set = features.omnians(network)

**Cut Edges and Vertices**

The :func:`phylozoo.core.network.sdnetwork.features.cut_edges` and :func:`phylozoo.core.network.sdnetwork.features.cut_vertices`
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
network types (galled, tree-child, tree-based), and structural properties
like binary resolution.

Basic Network Properties
^^^^^^^^^^^^^^^^

**Phylogenetic Trees**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_tree` function checks if the network
is a phylogenetic tree (no hybrid nodes).

.. code-block:: python

   from phylozoo.core.network.sdnetwork import classifications
   
   # Check if network is a tree
   is_tree = classifications.is_tree(network)

**Binary Networks**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_binary` function checks if the network
is binary (all internal nodes have degree 3, except for the root node which must have degree 2).

.. code-block:: python

   # Check if network is binary
   is_binary = classifications.is_binary(network)

**Simple Networks**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_simple` function checks if the network
is simple: containing exactly one non-leaf blob.

.. code-block:: python

   # Check if network is simple
   is_simple = classifications.is_simple(network)

**Parallel Edges**

The :func:`phylozoo.core.network.sdnetwork.classifications.has_parallel_edges` function checks if
the network contains parallel edges (multiple edges between the same pair of nodes).

.. code-block:: python

   # Check for parallel edges
   has_parallel = classifications.has_parallel_edges(network)

Network Classes
^^^^^^^^^^^^^^^

**Reticulation Number**

The :func:`phylozoo.core.network.sdnetwork.classifications.reticulation_number` function calculates
the reticulation number, which is the total number of hybrid edges minus the total number
of hybrid nodes.

.. code-block:: python

   # Calculate reticulation number
   reticulation_number = classifications.reticulation_number(network)

**Level and Vertex Level**

The :func:`phylozoo.core.network.sdnetwork.classifications.level` and :func:`phylozoo.core.network.sdnetwork.classifications.vertex_level` functions calculate the level and vertex level of the network, respectively.

The level is the maximum over all blobs of (number of hybrid edges minus number of hybrid nodes) in that blob.
The vertex level is the maximum over all blobs of the number of hybrid nodes in that blob.
Note that the vertex level is always less than or equal to the level and they coincide for binary networks.

.. code-block:: python

   # Calculate level and vertex level
   level = classifications.level(network)
   vertex_level = classifications.vertex_level(network)

**Stack-Free Networks**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_stackfree` function checks if the network
is stack-free, meaning no two hybrid nodes share a common parent.

.. code-block:: python

   # Check if network is stack-free
   is_stackfree = classifications.is_stackfree(network)

**Tree-Child Networks**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_strongly_treechild` function checks if the
network is strongly tree-child, meaning each internal node has at least one tree child. The
:func:`phylozoo.core.network.sdnetwork.classifications.is_weakly_treechild` function checks if the
network is weakly tree-child.

.. code-block:: python

   # Check if network is tree-child
   is_strongly_treechild = classifications.is_strongly_treechild(network)
   is_weakly_treechild = classifications.is_weakly_treechild(network)

**Tree-Based Networks**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_strongly_treebased` and :func:`phylozoo.core.network.sdnetwork.classifications.is_weakly_treebased` functions check if the
network is strongly tree-based and weakly tree-based, respectively.

.. code-block:: python

   # Check if network is tree-based
   is_strongly_treebased = classifications.is_strongly_treebased(network)
   is_weakly_treebased = classifications.is_weakly_treebased(network)

**Galled Networks**

The :func:`phylozoo.core.network.sdnetwork.classifications.is_galled` function checks if the network
is galled, meaning each reticulation is in its own cycle. Galled networks have special
properties that make them easier to analyze.

.. code-block:: python

   # Check if network is galled
   is_galled = classifications.is_galled(network)

Transformations
------------------------

The transformations module provides functions to modify network structures. These include
identifying parallel edges and suppressing 2-blobs.

**Parallel Edge Identification**

The :func:`phylozoo.core.network.sdnetwork.transformations.identify_parallel_edges` function
identifies parallel edges by replacing all parallel edges between the same pair of nodes with a single edge, and suppressing all resulting degree-2 nodes.

.. code-block:: python

   from phylozoo.core.network.sdnetwork import transformations
   
   # Identify parallel edges
   identified_net = transformations.identify_parallel_edges(network)

**Suppressing 2-Blobs**

The :func:`phylozoo.core.network.sdnetwork.transformations.suppress_2_blobs` function suppresses
all 2-blobs and the resulting degree-2 nodes.

.. code-block:: python

   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

Derivations
-------------------

The derivations module provides functions to extract derived structures from networks,
including conversion to directed networks, tree-of-blobs, subnetworks for specific
taxa sets, displayed trees, and induced splits and quartets.

Phylogenetic Trees
^^^^^^^^^^^^^^^^^^^

**Tree-of-Blobs**

The :func:`phylozoo.core.network.sdnetwork.derivations.tree_of_blobs` function extracts the
tree structure of blobs, representing the high-level topology of the network. This
simplifies the network by collapsing each blob into a single node.

.. code-block:: python

   from phylozoo.core.network.sdnetwork import derivations
   
   # Extract tree of blobs
   tob = derivations.tree_of_blobs(network)

**Displayed Trees**

The :func:`phylozoo.core.network.sdnetwork.derivations.displayed_trees` function generates all
displayed trees (trees embedded in the network). Optionally the function can return the probability of each displayed tree, calculated as the product of the gamma values of the hybrid edges contributing to the displayed tree.

.. code-block:: python

   # Get displayed trees
   for tree in derivations.displayed_trees(network):
       # Process each displayed tree
       pass

Phylogenetic Networks
^^^^^^^^^^^^^^^^^^^^^^

**Subnetworks**

The :func:`phylozoo.core.network.sdnetwork.derivations.subnetwork` function extracts a subnetwork
induced by a specific set of taxa. The subnetwork is defined as the union of all paths
from the requested leaves. The induced subgraph is taken, then degree-2 internal nodes are suppressed.

The function has three post-processing options: it can suppress 2-blobs, identify/merge parallel edges, and/or convert the result to a directed network.

.. code-block:: python

   # Get basic subnetwork
   subnetwork = derivations.subnetwork(network, taxa=["A", "B", "C"])

**k-Taxon Subnetworks**

The :func:`phylozoo.core.network.sdnetwork.derivations.k_taxon_subnetworks` function generates all
subnetworks induced by exactly k taxa. For each combination of k taxa, the corresponding
subnetwork is computed using the `subnetwork` function.

.. code-block:: python

   # Get all 4-taxon subnetworks
   k_subnets = derivations.k_taxon_subnetworks(network, k=4)

**Conversion to Directed Network**

The :func:`phylozoo.core.network.sdnetwork.derivations.to_d_network` function converts a
semi-directed network to a directed network by orienting all undirected edges away from a root.

.. code-block:: python

   # Convert to directed network
   d_net = derivations.to_d_network(network)

The :func:`phylozoo.core.network.sdnetwork.derivations.root_at_outgroup` function roots a
semi-directed network at the edge leading to a specified outgroup taxon (convenience wrapper
around :func:`to_d_network`).

.. code-block:: python

   # Root at outgroup leaf
   d_net = derivations.root_at_outgroup(network, outgroup='A')

Pairwise Distances
^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.network.sdnetwork.derivations.distances` function computes pairwise distances
between taxa. The function can compute shortest, longest, or
probability-weighted average distances across all displayed trees of the network.

.. code-block:: python

   # Compute distance matrix with average distances
   from phylozoo.core.network.sdnetwork.derivations import distances
   distance_matrix = distances(network, mode='average')  # or 'shortest', 'longest'

Partitions and Splits
^^^^^^^^^^^^^^^

**Partition from Blob**

The :func:`phylozoo.core.network.sdnetwork.derivations.partition_from_blob` function extracts the partition
of taxa induced by removing a specific blob from the network. 

.. code-block:: python

   # Extract partition from specific blob
   from phylozoo.core.network.sdnetwork.features import blobs
   blob = list(blobs(network))[0]
   partition = derivations.partition_from_blob(network, blob)

**Split from Cut-Edge**

The :func:`phylozoo.core.network.sdnetwork.derivations.split_from_cutedge` function extracts the split
induced by a specific cut edge. The split is the 2-partition of taxa obtained when removing the edge from the network.

.. code-block:: python

   # Extract split from specific cut edge
   from phylozoo.core.network.sdnetwork.features import cut_edges
   cut_edge = list(cut_edges(network))[0]
   split = derivations.split_from_cutedge(network, cut_edge)

**Induced Splits**

The :func:`phylozoo.core.network.sdnetwork.derivations.induced_splits` function extracts all splits
induced by cut-edges of the network. This returns a :class:`phylozoo.core.split.SplitSystem` object.

.. code-block:: python

   # Extract induced splits
   splits = derivations.induced_splits(network)

**Displayed Splits**

The :func:`phylozoo.core.network.sdnetwork.derivations.displayed_splits` function extracts all splits
induced by all displayed trees, weighted by their probabilities. This returns a :class:`phylozoo.core.split.WeightedSplitSystem` object.

.. code-block:: python

   # Extract displayed splits
   displayed_splits = derivations.displayed_splits(network)

Quartets
^^^^^^^^

The :func:`phylozoo.core.network.sdnetwork.derivations.displayed_quartets` function extracts quartet
profiles from displayed trees. This returns a :class:`phylozoo.core.quartet.QuartetProfileSet` object.

.. code-block:: python

   # Extract displayed quartets
   quartets = derivations.displayed_quartets(network)


Isomorphism Checking
--------------------

The isomorphism module provides functions to check if two networks have the same
topological structure. 

The :func:`phylozoo.core.network.sdnetwork.isomorphism.is_isomorphic` function checks if two
networks are isomorphic. Labels are always checked, and additional attributes can
be specified for comparison.

.. code-block:: python

   from phylozoo.core.network.sdnetwork import isomorphism
   
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

- :doc:`API Reference <../../../api/core/network>` - Complete function signatures and detailed examples
- :doc:`Semi-Directed Network Class <semi_directed_network_class>` - The SemiDirectedPhyNetwork class and its properties
- :doc:`Semi-Directed Generator <generators>` - Level-k network generators
- :doc:`Directed Networks <../directed/overview>` - Fully directed network representations
