Directed Networks (Advanced)
=============================

The :mod:`phylozoo.core.network.dnetwork` module provides advanced features, transformations,
and analysis capabilities for :class:`DirectedPhyNetwork`. This page covers network features,
classifications, transformations, derivations, and isomorphism checking. For basic network
operations, see :doc:`Directed Networks (Basic) <basic>`.

All classes and functions on this page can be imported from the core network module:

.. code-block:: python

   from phylozoo.core.dnetwork import features, classifications, transformations, derivations, isomorphism

Network Features
----------------

The features module provides functions to extract advanced structural properties of
directed networks. These include LSA (Least Stable Ancestor) nodes, blobs (maximal
biconnected components), omnian nodes, and cut edges/vertices that characterize
network topology.

Structural Analysis
^^^^^^^^^^^^^^^^^^

**LSA Node**

The :func:`phylozoo.core.dnetwork.features.lsa_node` function finds the Least Stable
Ancestor node, which is the lowest node through which all root-to-leaf paths pass.
The LSA is a fundamental concept in phylogenetic network analysis, representing the
most recent common ancestor of all leaves.

.. code-block:: python

   from phylozoo.core.dnetwork import features
   
   # Find LSA (Least Stable Ancestor) node
   lsa = features.lsa_node(network)

**Blobs**

The :func:`phylozoo.core.dnetwork.features.blobs` function extracts all blobs (maximal
biconnected components) in the network. Blobs represent the reticulate structure of
the network, where cycles and hybrid nodes create biconnected regions.

.. code-block:: python

   # Get blobs (maximal biconnected components)
   blob_list = features.blobs(network, trivial=False, leaves=False)

**K-Blobs**

The :func:`phylozoo.core.dnetwork.features.k_blobs` function gets blobs with exactly
k incident edges, useful for identifying cycles and reticulation structures. For example,
2-blobs are cycles in the network.

.. code-block:: python

   # Get k-blobs (blobs with exactly k incident edges)
   k2_blobs = features.k_blobs(network, k=2)

**Omnians**

The :func:`phylozoo.core.dnetwork.features.omnians` function finds omnian nodes,
which are nodes that are ancestors of all leaves. These nodes are important for
understanding network structure and root placement.

.. code-block:: python

   # Find omnians (nodes that are ancestors of all leaves)
   omnians_set = features.omnians(network)

**Cut Edges and Vertices**

The :func:`phylozoo.core.dnetwork.features.cut_edges` and :func:`phylozoo.core.dnetwork.features.cut_vertices`
functions identify edges and vertices whose removal disconnects the network. These
are important for understanding network connectivity and decomposition.

.. code-block:: python

   # Get cut edges and cut vertices
   cut_edges = features.cut_edges(network)
   cut_vertices = features.cut_vertices(network)

Network Classifications
-----------------------

The classifications module provides functions to determine various mathematical and
phylogenetic properties of networks. These include level (number of reticulations),
network types (LSA, galled, tree-child, tree-based, normal), and structural properties
like binary resolution and ultrametricity.

Basic Properties
^^^^^^^^^^^^^^^^

**Level**

The :func:`phylozoo.core.dnetwork.classifications.level` function calculates the level
(number of reticulations) of the network. The level is a fundamental measure of network
complexity, representing the minimum number of edges that must be removed to obtain a tree.

.. code-block:: python

   from phylozoo.core.dnetwork import classifications
   
   # Calculate level (number of reticulations)
   level = classifications.level(network)

**Binary Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_binary` function checks if the network
is binary (all internal nodes have degree 3). Binary networks are a common assumption
in many phylogenetic algorithms.

.. code-block:: python

   # Check if network is binary
   is_binary = classifications.is_binary(network)

**Ultrametric Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_ultrametric` function checks if all
root-to-leaf distances are equal, which is important for certain types of phylogenetic
analyses.

.. code-block:: python

   # Check if network is ultrametric
   is_ultrametric = classifications.is_ultrametric(network)

Network Types
^^^^^^^^^^^^^

**LSA Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_lsa_network` function checks if the
network is an LSA network, meaning the root equals the LSA node. This is a desirable
property for many network analyses.

.. code-block:: python

   # Check if network is LSA network
   is_lsa = classifications.is_lsa_network(network)

**Galled Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_galled` function checks if the network
is galled, meaning each reticulation is in its own cycle. Galled networks have special
properties that make them easier to analyze.

.. code-block:: python

   # Check if network is galled
   is_galled = classifications.is_galled(network)

**Tree-Child Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_treechild` function checks if the
network is tree-child, meaning each internal node has at least one tree child. This
property ensures that the network has a tree-like structure.

.. code-block:: python

   # Check if network is tree-child
   is_treechild = classifications.is_treechild(network)

**Tree-Based Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_treebased` function checks if the
network is tree-based, meaning it can be obtained from a tree by adding edges. Tree-based
networks are a broad class that includes many biologically relevant networks.

.. code-block:: python

   # Check if network is tree-based
   is_treebased = classifications.is_treebased(network)

**Normal Networks**

The :func:`phylozoo.core.dnetwork.classifications.is_normal` function checks if the network
is normal, meaning it is tree-child without shortcuts. Normal networks have additional
structural constraints.

.. code-block:: python

   # Check if network is normal
   is_normal = classifications.is_normal(network)

**Parallel Edges**

The :func:`phylozoo.core.dnetwork.classifications.has_parallel_edges` function checks if
the network contains parallel edges (multiple edges between the same pair of nodes).

.. code-block:: python

   # Check for parallel edges
   has_parallel = classifications.has_parallel_edges(network)

Advanced Transformations
------------------------

The transformations module provides functions to modify network structures while
preserving essential topological information. These include converting to LSA form,
binary resolution of high-degree nodes, identifying parallel edges, and suppressing
degree-2 nodes in 2-blobs.

**LSA Network Conversion**

The :func:`phylozoo.core.dnetwork.transformations.to_lsa_network` function converts
a network to LSA form by suppressing nodes above the LSA. This transformation simplifies
the network while preserving its essential structure.

.. code-block:: python

   from phylozoo.core.dnetwork import transformations
   
   # Convert to LSA network
   lsa_net = transformations.to_lsa_network(network)

**Binary Resolution**

The :func:`phylozoo.core.dnetwork.transformations.binary_resolution` function resolves
high-degree nodes to binary form using caterpillar structures, preserving gamma values
and branch lengths. This transformation is useful for algorithms that require binary networks.

.. code-block:: python

   # Binary resolution (resolve high-degree nodes)
   binary_net = transformations.binary_resolution(network)

**Parallel Edge Identification**

The :func:`phylozoo.core.dnetwork.transformations.identify_parallel_edges` function
identifies parallel edges by adding explicit edge keys. This transformation makes
parallel edges explicit in the network representation.

.. code-block:: python

   # Identify parallel edges
   identified_net = transformations.identify_parallel_edges(network)

**Suppressing 2-Blobs**

The :func:`phylozoo.core.dnetwork.transformations.suppress_2_blobs` function suppresses
degree-2 nodes in 2-blobs, simplifying network structure while preserving topology.

.. code-block:: python

   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

Network Derivations
-------------------

The derivations module provides functions to extract derived structures from networks,
including conversion to semi-directed networks, tree of blobs, subnetworks for specific
taxa sets, displayed trees, and induced splits and quartets. These operations are
fundamental for network-based phylogenetic analysis.

**Conversion to Semi-Directed Network**

The :func:`phylozoo.core.dnetwork.derivations.to_sd_network` function converts a
directed network to a semi-directed network by undirecting non-hybrid edges. This
allows you to work with semi-directed network algorithms on directed networks.

.. code-block:: python

   from phylozoo.core.dnetwork import derivations
   
   # Convert to semi-directed network
   sd_net = derivations.to_sd_network(network)

**Tree of Blobs**

The :func:`phylozoo.core.dnetwork.derivations.tree_of_blobs` function extracts the
tree structure of blobs, representing the high-level topology of the network. This
simplifies the network by collapsing each blob into a single node.

.. code-block:: python

   # Extract tree of blobs
   blob_tree = derivations.tree_of_blobs(network)

**Subnetworks**

The :func:`phylozoo.core.dnetwork.derivations.subnetwork` function extracts a subnetwork
induced by a specific set of taxa, and :func:`phylozoo.core.dnetwork.derivations.k_taxon_subnetworks`
generates all k-taxon subnetworks. These operations are useful for analyzing relationships
among specific subsets of taxa.

.. code-block:: python

   # Get subnetwork for specific taxa
   subnetwork = derivations.subnetwork(network, taxa={"A", "B", "C"})
   
   # Get k-taxon subnetworks
   k_subnets = derivations.k_taxon_subnetworks(network, k=4)

**Displayed Trees**

The :func:`phylozoo.core.dnetwork.derivations.displayed_trees` function generates all
displayed trees (trees embedded in the network), optionally including edge probabilities.
Displayed trees represent all possible tree topologies that can be extracted from the network.

.. code-block:: python

   # Get displayed trees
   for tree in derivations.displayed_trees(network):
       # Process each displayed tree
       pass

**Splits and Quartets**

The :func:`phylozoo.core.dnetwork.derivations.induced_splits` function extracts splits
induced by the network, and :func:`phylozoo.core.dnetwork.derivations.displayed_quartets`
extracts quartet profiles from displayed trees. These operations connect network structures
to split and quartet-based phylogenetic methods.

.. code-block:: python

   # Extract splits and quartets
   splits = derivations.induced_splits(network)
   quartets = derivations.displayed_quartets(network)

Isomorphism Checking
--------------------

The isomorphism module provides functions to check if two networks have the same
topological structure. Isomorphism checking considers node labels by default and can
optionally compare additional node, edge, and graph attributes to determine if networks
are structurally equivalent.

**Isomorphism**

The :func:`phylozoo.core.dnetwork.isomorphism.is_isomorphic` function checks if two
networks are isomorphic. Labels are always checked, and additional attributes can
be specified for comparison. This is useful for identifying duplicate networks or
verifying network transformations.

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

.. note::
   All transformation and derivation functions return new network instances. Networks 
   are immutable.

.. tip::
   Use network classifications to determine which algorithms can be applied to your 
   network. For example, certain diversity optimization methods require networks 
   without 2-blobs or parallel edges.

.. warning::
   Some advanced operations (like binary resolution) can significantly increase network 
   size. Use with caution on large networks.

See Also
--------

- :doc:`Directed Networks (Basic) <basic>` - Basic network operations
- :doc:`Semi-Directed Networks <../semi_directed/overview>` - Semi-directed network representations
