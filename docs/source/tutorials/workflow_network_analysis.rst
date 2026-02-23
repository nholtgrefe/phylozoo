Network Analysis Workflow
==========================

This page demonstrates a complete workflow for analyzing phylogenetic networks using 
PhyloZoo. We'll use a consistent example network throughout to show how different 
operations work together.

Example Network
---------------

We'll use a simple 4-taxon network with one hybrid node:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Create example network
   network = DirectedPhyNetwork(
       edges=[
           ("root", "u1"), ("root", "u2"),  # Root to tree nodes
           ("u1", "h", {"gamma": 0.6}),     # Hybrid edge 1
           ("u2", "h", {"gamma": 0.4}),     # Hybrid edge 2 (sum = 1.0)
           ("h", "A"), ("u1", "B"), ("u2", "C"), ("root", "D")
       ],
       nodes=[
           ("A", {"label": "A"}),
           ("B", {"label": "B"}),
           ("C", {"label": "C"}),
           ("D", {"label": "D"})
       ]
   )

Step 1: Load or Create Network
-------------------------------

Networks can be created programmatically (as above) or loaded from files:

.. code-block:: python

   # Load from file
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Or create from string
   network = DirectedPhyNetwork.from_string(
       "(A:0.5, B:0.3, (C:0.2, D:0.4)h:0.1);",
       format="enewick"
   )

See :doc:`I/O <../manual/io>` for supported formats.

Step 2: Basic Network Properties
--------------------------------

Examine basic network structure:

.. code-block:: python

   # Basic counts
   print(f"Nodes: {network.num_nodes}, Edges: {network.num_edges}")
   print(f"Root: {network.root_node}")
   print(f"Leaves: {network.leaves}")
   print(f"Taxa: {network.taxa}")
   
   # Node types
   print(f"Hybrid nodes: {network.hybrid_nodes}")
   print(f"Tree nodes: {network.tree_nodes}")
   
   # Check if tree
   print(f"Is tree: {network.is_tree()}")

Step 3: Network Classifications
-------------------------------

Classify the network:

.. code-block:: python

   from phylozoo.core.dnetwork import classifications
   
   # Basic properties
   level = classifications.level(network)
   is_binary = classifications.is_binary(network)
   is_lsa = classifications.is_lsa_network(network)
   
   print(f"Level: {level}")
   print(f"Binary: {is_binary}")
   print(f"LSA network: {is_lsa}")

Step 4: Extract Network Features
--------------------------------

Extract advanced features:

.. code-block:: python

   from phylozoo.core.dnetwork import features
   
   # Find LSA node
   lsa = features.lsa_node(network)
   
   # Get blobs
   blobs = features.blobs(network, trivial=False, leaves=False)
   
   # Get k-blobs
   k2_blobs = features.k_blobs(network, k=2)
   
   print(f"LSA node: {lsa}")
   print(f"Number of blobs: {len(blobs)}")

Step 5: Extract Derived Structures
-----------------------------------

Extract splits, quartets, and trees from the network:

.. code-block:: python

   from phylozoo.core.dnetwork import derivations
   
   # Extract splits
   splits = derivations.induced_splits(network)
   print(f"Number of splits: {len(splits)}")
   
   # Extract quartets
   quartets = derivations.displayed_quartets(network)
   print(f"Number of quartet profiles: {len(quartets)}")
   
   # Get displayed trees
   trees = list(derivations.displayed_trees(network))
   print(f"Number of displayed trees: {len(trees)}")
   
   # Extract subnetwork
   subnetwork = derivations.subnetwork(network, taxa={"A", "B", "C"})

Step 6: Network Transformations
-------------------------------

Transform the network:

.. code-block:: python

   from phylozoo.core.dnetwork import transformations
   
   # Convert to LSA network
   lsa_net = transformations.to_lsa_network(network)
   
   # Binary resolution
   binary_net = transformations.binary_resolution(network)
   
   # Suppress 2-blobs
   simplified = transformations.suppress_2_blobs(network)

Step 7: Convert to Semi-Directed
--------------------------------

Convert to semi-directed representation:

.. code-block:: python

   from phylozoo.core.dnetwork.derivations import to_sd_network
   
   # Convert to semi-directed
   sd_network = to_sd_network(network)
   
   # Now can use semi-directed features
   from phylozoo.core.sdnetwork import features as sd_features
   blobs = sd_features.blobs(sd_network)
   root_locs = sd_features.root_locations(sd_network)

Step 8: Visualize Network
-------------------------

Visualize the network:

.. code-block:: python

   from phylozoo.viz.dnetwork import plot_network
   
   # Basic plot
   plot_network(network, show=True)
   
   # Save to file
   fig = plot_network(network, show=False)
   fig.savefig("network.png", dpi=300)

See :doc:`Visualization <../manual/viz>` for advanced plotting options.

Step 9: Save Results
--------------------

Save the network and extracted data:

.. code-block:: python

   # Save network
   network.save("analyzed_network.enewick")
   
   # Save splits
   splits.save("splits.nexus")
   
   # Save quartets (if supported)
   # quartets.save("quartets.txt")

Complete Example
----------------

Here's the complete workflow in one script:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.core.dnetwork import (
       classifications, features, transformations, derivations
   )
   from phylozoo.viz.dnetwork import plot_network
   
   # Load network
   network = DirectedPhyNetwork.load("input.enewick")
   
   # Analyze
   print(f"Level: {classifications.level(network)}")
   print(f"LSA: {features.lsa_node(network)}")
   
   # Extract structures
   splits = derivations.induced_splits(network)
   quartets = derivations.displayed_quartets(network)
   
   # Transform
   binary_net = transformations.binary_resolution(network)
   
   # Visualize
   plot_network(network, show=True)
   
   # Save
   network.save("output.enewick")
   splits.save("splits.nexus")

.. note::
   This workflow demonstrates basic to intermediate operations. For advanced features, 
   see :doc:`Networks (Advanced) <../manual/core/networks/advanced>`.

.. tip::
   Use network classifications to determine which algorithms can be applied. For 
   example, certain diversity optimization methods require networks without 2-blobs.

.. seealso::
   For detailed function documentation, see :doc:`Networks (Basic) <../manual/core/networks/basic>` 
   and :doc:`Networks (Advanced) <../manual/core/networks/advanced>`.
