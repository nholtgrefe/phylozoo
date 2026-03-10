Quickstart
==========

Get started with PhyloZoo in minutes! This guide walks you through installation,
creating networks, and working with related data types, at every step referring to the corresponding section in the manual.

Getting Started
---------------

Install PhyloZoo using pip. For most users (including visualization), use:

.. code-block:: bash

   pip install phylozoo[viz]

For a minimal install without plotting, use ``pip install phylozoo``.
See the :doc:`Installation Guide <../manual/installation>` for full details.

Networks
--------

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

This creates a simple tree with three leaves.

Now let's create a network with a hybrid node. Hybrid edges use a ``gamma`` value
(must sum to 1.0 for edges into the same hybrid):

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

PhyloZoo also supports **semi-directed networks**, which allow undirected tree edges
for modeling root uncertainty. See :doc:`Semi-Directed Networks <../manual/core/networks/semi_directed/overview>` for details.

Load and Save
~~~~~~~~~~~~~

Networks can be saved and loaded in standard formats:

.. code-block:: python

   # Save network to file (eNewick format)
   network.save("my_network.enewick")

   # Load network from file
   loaded_net = DirectedPhyNetwork.load("my_network.enewick")

Supported formats include eNewick, NEXUS, and DOT. See the :doc:`I/O overview <../manual/utils/io/index>` for all formats.

Visualization
~~~~~~~~~~~~~

Plot networks with the built-in visualization module:

.. code-block:: python

   from phylozoo.viz import plot

   # Plot network (opens interactive window)
   plot(network, show=True)

   # Save to file instead
   plot(network, path="network.png")

For layout options and styling, see :doc:`Visualization <../manual/visualization/overview>`.

Related Data Types
------------------

PhyloZoo provides integrated support for related data types used in phylogenetic
analysis.

Distance Matrices
~~~~~~~~~~~~~~~~~

Create distance matrices from NumPy arrays or compute them from sequences:

.. code-block:: python

   import numpy as np
   from phylozoo.core.distance import DistanceMatrix

   # Create from array
   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=["A", "B", "C"])

Sequences
~~~~~~~~~

Create multiple sequence alignments (MSAs) and compute Hamming distances:

.. code-block:: python

   from phylozoo import MSA
   from phylozoo.core.sequence import hamming_distances

   sequences = {"taxon1": "ACGTACGT", "taxon2": "TGCAACGT", "taxon3": "ACGTTTAA"}
   msa = MSA(sequences)
   dm = hamming_distances(msa)

See :doc:`Distance Matrices <../manual/core/distance>` and :doc:`Sequences <../manual/core/sequences>` for more.

Splits and Quartets
~~~~~~~~~~~~~~~~~~~

For split-based and quartet-based analyses:

.. code-block:: python

   from phylozoo.core.split import Split, SplitSystem
   from phylozoo.core.quartet import Quartet, QuartetProfile

   # Create a split (bipartition of taxa)
   split = Split({"A", "B"}, {"C", "D"})
   system = SplitSystem([split])

   # Create a quartet
   quartet = Quartet(split)
   profile = QuartetProfile([quartet])

Splits represent bipartitions of taxa; quartets are four-taxon unrooted trees.
See :doc:`Splits <../manual/core/splits/overview>` and :doc:`Quartets <../manual/core/quartets/overview>` for details.

Quick Links
-----------

* :doc:`Installation Guide <../manual/installation>` — Detailed installation and requirements
* :doc:`Introduction <../manual/introduction>` — PhyloZoo's design and features
* :doc:`Networks <../manual/core/networks/overview>` — Directed and semi-directed networks
* :doc:`Distance Matrices <../manual/core/distance>` — Pairwise distances
* :doc:`Sequences <../manual/core/sequences>` — Multiple sequence alignments
* :doc:`I/O <../manual/utils/io/index>` — File formats and operations
* :doc:`Visualization <../manual/visualization/overview>` — Plotting and styling

Next Steps
----------

Ready to dive deeper? Check out:

1. The :doc:`Manual <../manual/index>` for comprehensive guides on all modules
2. The :doc:`API Reference <../api/index>` for complete function signatures and examples
