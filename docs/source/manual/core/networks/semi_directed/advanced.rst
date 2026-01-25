Semi-Directed Networks (Advanced)
===================================

The :mod:`phylozoo.core.network.sdnetwork` module provides advanced features, transformations,
and analysis capabilities for :class:`SemiDirectedPhyNetwork`. This page covers network features,
classifications, transformations, and derivations. For basic network operations, see
:doc:`Semi-Directed Networks (Basic) <basic>`.

All classes and functions on this page can be imported from the core network module:

.. code-block:: python

   from phylozoo.core.sdnetwork import features, classifications, transformations, derivations

Network Features
----------------

The features module provides functions to extract advanced structural properties of
semi-directed networks. These include blobs (maximal biconnected components), root
locations (where the network can be rooted), and up-down paths that characterize
network topology.

**Blobs**

The :func:`phylozoo.core.sdnetwork.features.blobs` function extracts all blobs (maximal
biconnected components) in the network. Blobs represent the reticulate structure of
the network.

**Root Locations**

The :func:`phylozoo.core.sdnetwork.features.root_locations` function finds all possible
root locations (nodes or edges) where the network can be rooted.

**Up-Down Paths**

The :func:`phylozoo.core.sdnetwork.m_multigraph.features.updown_path_vertices` function
finds up-down paths between nodes, which are paths that go up (against directed edges)
and then down (with directed edges or along undirected edges).

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
-----------------------

The classifications module provides functions to determine various mathematical and
phylogenetic properties of semi-directed networks. These include level (number of
reticulations), network types (galled), and structural properties like binary resolution.

**Level**

The :func:`phylozoo.core.sdnetwork.classifications.level` function calculates the level
(number of reticulations) of the network.

**Network Types**

The :func:`phylozoo.core.sdnetwork.classifications.is_galled` function checks if the
network is galled (no hybrid node is ancestral to another hybrid node in the same blob).

**Structural Properties**

The :func:`phylozoo.core.sdnetwork.classifications.is_binary` function checks if the network
is binary (all internal nodes have degree 3).

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

The transformations module provides functions to modify network structures while
preserving essential topological information.

**Suppressing 2-Blobs**

The :func:`phylozoo.core.sdnetwork.transformations.suppress_2_blobs` function removes
degree-2 nodes in 2-blobs, simplifying the network structure without changing its
essential topology.

.. code-block:: python

   from phylozoo.core.sdnetwork import transformations
   
   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

Network Derivations
-------------------

The derivations module provides functions to extract derived structures from semi-directed
networks, including tree of blobs and subnetworks for specific taxa sets. These operations
are fundamental for network-based phylogenetic analysis and allow you to focus on specific
subsets of taxa or extract simplified representations of network structure.

**Tree of Blobs**

The :func:`phylozoo.core.sdnetwork.derivations.tree_of_blobs` function extracts the
tree structure of blobs, representing the high-level topology of the network.

**Subnetworks**

The :func:`phylozoo.core.sdnetwork.derivations.subnetwork` function extracts a subnetwork
induced by a specific set of taxa, and :func:`phylozoo.core.sdnetwork.derivations.k_taxon_subnetworks`
generates all k-taxon subnetworks.

.. code-block:: python

   from phylozoo.core.sdnetwork import derivations
   
   # Extract tree of blobs
   blob_tree = derivations.tree_of_blobs(network)
   
   # Get subnetwork for specific taxa
   subnetwork = derivations.subnetwork(network, taxa={"A", "B", "C"})
   
   # Get k-taxon subnetworks
   k_subnets = derivations.k_taxon_subnetworks(network, k=4)

.. note::
   All transformation and derivation functions return new network instances. Networks 
   are immutable.

.. tip::
   Use network classifications to determine which algorithms can be applied to your 
   network.

See Also
--------

- :doc:`Semi-Directed Networks (Basic) <basic>` - Basic network operations
- :doc:`Directed Networks <../directed/overview>` - Fully directed network representations
