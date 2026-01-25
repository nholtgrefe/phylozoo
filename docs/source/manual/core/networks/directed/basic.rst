Directed Networks (Basic)
==========================

The :mod:`phylozoo.core.network.dnetwork` module provides the :class:`DirectedPhyNetwork` class,
which represents fully directed phylogenetic networks with a single root node. DirectedPhyNetwork
is fundamental to phylogenetic analysis, explicitly representing the direction of evolutionary
time through fully directed edges. All edges are directed, and hybrid nodes have in-degree >= 2
and out-degree 1. For advanced features like network analysis, transformations, and classifications,
see :doc:`Directed Networks (Advanced) <advanced>`.

All classes and functions on this page can be imported from the core network module:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   # or directly
   from phylozoo.core.network.dnetwork import DirectedPhyNetwork

Working with DirectedPhyNetwork
--------------------------------

The :class:`phylozoo.core.network.dnetwork.DirectedPhyNetwork` class is the canonical container
for fully directed phylogenetic networks in PhyloZoo. It provides an immutable representation
that ensures data integrity throughout your analysis pipeline.

.. note::
   :class: dropdown

   **Implementation details**

   DirectedPhyNetwork is designed for immutability and performance:

   - Networks are immutable after construction
   - All edges are directed, explicitly representing time flow
   - Hybrid nodes are validated to have in-degree >= 2 and out-degree 1
   - Edge attributes (branch lengths, bootstrap, gamma) are stored efficiently
   - The network structure is validated at construction time

   For implementation details, see :mod:`src/phylozoo/core/network/dnetwork/base.py`.

.. figure:: ../../../../_static/images/example_tree_directed.png
   :alt: Example directed tree
   :width: 400px
   :align: center
   
   Example of a simple directed tree.

.. figure:: ../../../../_static/images/example_hybrid_directed.png
   :alt: Example directed network with hybrid
   :width: 400px
   :align: center
   
   Example of a directed network with a hybrid node.

Creating a DirectedPhyNetwork
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DirectedPhyNetwork objects can be created from edge and node specifications. The constructor
accepts a list of edges (as tuples) and optional node attributes. Each edge can have
attributes like branch lengths, bootstrap values, and gamma probabilities for hybrid edges.

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Simple tree
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B"), ("root", "C")],
       nodes=[
           ("A", {"label": "A"}),
           ("B", {"label": "B"}),
           ("C", {"label": "C"})
       ]
   )

Networks with hybrid nodes require multiple incoming edges to the hybrid node, each with
a gamma probability. The gamma values for all edges entering a hybrid node must sum to 1.0.

.. code-block:: python

   # Network with hybridization
   hybrid_net = DirectedPhyNetwork(
       edges=[
           ("root", "u1"), ("root", "u2"),
           ("u1", "h", {"gamma": 0.6}),
           ("u2", "h", {"gamma": 0.4}),  # Must sum to 1.0
           ("h", "leaf1")
       ],
       nodes=[("leaf1", {"label": "A"})]
   )

Accessing Network Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DirectedPhyNetwork provides properties to access fundamental network structure information,
including node and edge counts, root and leaf identification, and node type classification.

.. code-block:: python

   # Node and edge counts
   num_nodes = network.num_nodes
   num_edges = network.num_edges
   
   # Network structure
   root = network.root_node
   leaves = network.leaves  # Set of leaf node IDs
   taxa = network.taxa       # Set of taxon labels (strings)
   internal_nodes = network.internal_nodes
   
   # Node types
   hybrid_nodes = network.hybrid_nodes
   tree_nodes = network.tree_nodes
   
   # Check if tree
   is_tree = network.is_tree()

Accessing Edge Attributes
^^^^^^^^^^^^^^^^^^^^^^^^^

Edge attributes such as branch lengths, bootstrap support values, and gamma probabilities
can be accessed using dedicated methods. Branch lengths represent evolutionary distances,
bootstrap values indicate statistical support, and gamma probabilities specify the
contribution of each parent edge to a hybrid node.

.. code-block:: python

   # Get branch length
   bl = network.get_branch_length("root", "A")
   
   # Get bootstrap support
   bs = network.get_bootstrap("root", "A")
   
   # Get gamma (for hybrid edges)
   gamma = network.get_gamma("u1", "h")

File Input/Output
^^^^^^^^^^^^^^^^^

DirectedPhyNetwork supports loading from and saving to files in multiple formats. The
default format is eNewick, which extends the standard Newick format to support hybrid
nodes and edge attributes.

**Supported Formats:**

- **eNewick** (default): Extended Newick format. Extensions: ``.enewick``, ``.eNewick``, ``.enwk``, ``.nwk``, ``.newick``
- **DOT**: Graphviz format. Extensions: ``.dot``, ``.gv``

.. code-block:: python

   # Save network
   network.save("network.enewick")
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")

.. seealso::
   The I/O system uses the :class:`phylozoo.utils.io.IOMixin` interface, providing
   consistent file handling across PhyloZoo classes. For details on the I/O system,
   see the :doc:`I/O documentation <../../../../io>`. For specific information about
   supported file formats and parameter options for networks, see the
   :mod:`API reference <phylozoo.core.network.dnetwork.io>`.

Basic Transformations
---------------------

Basic transformation functions allow you to simplify network structures while preserving
topological information.

**Suppressing 2-Blobs**

The :func:`phylozoo.core.dnetwork.transformations.suppress_2_blobs` function removes
degree-2 nodes within 2-blobs (biconnected components with 2 incident edges), simplifying
the network without changing its essential structure.

.. code-block:: python

   from phylozoo.core.dnetwork import transformations
   
   # Suppress degree-2 nodes (simplify network)
   simplified = transformations.suppress_2_blobs(network)

For more advanced transformations, see :doc:`Directed Networks (Advanced) <advanced>`.

.. note::
   For advanced network features like LSA nodes, blobs, network classifications, 
   binary resolution, and isomorphism checking, see :doc:`Directed Networks (Advanced) <advanced>`.

.. tip::
   All networks are immutable. To modify a network, create a new instance with the 
   desired changes.

See Also
--------

- :doc:`Directed Networks (Advanced) <advanced>` - Advanced features, transformations, and classifications
- :doc:`Semi-Directed Networks <../semi_directed/overview>` - Semi-directed network representations
- :doc:`I/O <../../../../io>` - File I/O operations and formats
