Quickstart
==========

Get started with PhyloZoo in minutes! This guide will help you install PhyloZoo and 
create your first phylogenetic network.

Installation
------------

Install PhyloZoo using pip:

.. code-block:: bash

   pip install phylozoo

For development installation:

.. code-block:: bash

   git clone https://github.com/yourusername/phylozoo.git
   cd phylozoo
   pip install -e ".[dev]"

Requirements: Python >= 3.10, NumPy >= 1.20.0, NetworkX >= 3.0.0

Your First Network
------------------

Let's create a simple phylogenetic network:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Create a simple tree network
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B"), ("root", "C")],
       nodes=[
           ("A", {"label": "A"}),
           ("B", {"label": "B"}),
           ("C", {"label": "C"})
       ]
   )
   
   print(f"Network has {network.num_nodes} nodes")
   print(f"Leaves: {network.leaves}")
   print(f"Is tree: {network.is_tree()}")

This creates a simple tree with three leaves. Now let's create a network with a hybrid node:

.. code-block:: python

   # Create a network with hybridization
   hybrid_net = DirectedPhyNetwork(
       edges=[
           ("root", "u1"), ("root", "u2"),  # Root to tree nodes
           ("u1", "h", {"gamma": 0.6}),    # Hybrid edge
           ("u2", "h", {"gamma": 0.4}),     # Hybrid edge (must sum to 1.0)
           ("h", "leaf1")                   # Hybrid to leaf
       ],
       nodes=[("leaf1", {"label": "A"})]
   )
   
   print(f"Network level: {hybrid_net.level()}")

Basic Operations
----------------

**Load and Save Networks**

.. code-block:: python

   # Save network to file
   network.save("my_network.enewick")
   
   # Load network from file
   loaded_net = DirectedPhyNetwork.load("my_network.enewick")

See the :doc:`I/O <manual/io>` page for supported formats.

**Visualize Networks**

.. code-block:: python

   from phylozoo.viz.dnetwork import plot_network
   
   # Plot network
   plot_network(network, show=True)

**Work with Sequences**

.. code-block:: python

   from phylozoo import MSA
   
   # Create MSA from sequences
   sequences = {
       "taxon1": "ACGTACGT",
       "taxon2": "ACGTACGT",
       "taxon3": "ACGTTTAA"
   }
   msa = MSA(sequences)
   
   # Compute distance matrix
   from phylozoo.core.sequence import hamming_distances
   distances = hamming_distances(msa)

Quick Links
-----------

* :doc:`Installation Guide <manual/installation>` - Detailed installation instructions
* :doc:`Introduction <manual/introduction>` - Learn about PhyloZoo's design and features
* :doc:`Networks (Basic) <manual/core/networks/basic>` - Basic network operations
* :doc:`Network Analysis Workflow <../tutorials/workflow_network_analysis>` - Complete workflow example
* :doc:`API Reference <library/index>` - Full API documentation

Next Steps
----------

Ready to dive deeper? Check out:

1. The :doc:`Manual <manual/index>` for comprehensive guides
2. :doc:`Tutorials <tutorials/index>` for extended examples
3. The :doc:`Library <library/index>` for complete API reference

For help and support, see the :doc:`Guides <guides/index>` section.
