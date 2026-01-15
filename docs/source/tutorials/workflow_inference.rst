Inference Workflow
==================

This page demonstrates a complete workflow for inferring phylogenetic networks from 
sequence data using PhyloZoo. We'll use consistent example data throughout to show 
how the inference pipeline works.

Example Data
------------

We'll use a simple 4-taxon MSA:

.. code-block:: python

   from phylozoo import MSA
   
   # Create example MSA
   sequences = {
       "A": "ACGTACGTACGT",
       "B": "ACGTACGTACGT",
       "C": "ACGTTTAAACGT",
       "D": "ACGTTTAAACGT"
   }
   msa = MSA(sequences)

Step 1: Load Sequence Data
--------------------------

Load sequences from files:

.. code-block:: python

   # Load from FASTA
   msa = MSA.load("alignment.fasta")
   
   # Or from NEXUS
   msa = MSA.load("alignment.nexus", format="nexus")

See :doc:`I/O <../manual/io>` for supported formats.

Step 2: Compute Distance Matrix
--------------------------------

Compute distance matrix from sequences:

.. code-block:: python

   from phylozoo.core.sequence import hamming_distances
   
   # Compute Hamming distances
   distances = hamming_distances(msa)
   
   print(f"Distance A-B: {distances.get_distance('A', 'B')}")
   print(f"Distance A-C: {distances.get_distance('A', 'C')}")

Alternatively, load a pre-computed distance matrix:

.. code-block:: python

   from phylozoo import DistanceMatrix
   distances = DistanceMatrix.load("distances.nexus")

Step 3: Generate Quartet Profiles
----------------------------------

Generate quartet profiles from distance matrix using delta heuristic:

.. code-block:: python

   from phylozoo.inference.squirrel import delta_heuristic
   
   # Generate quartet profiles
   profile_set = delta_heuristic(distances)
   
   print(f"Number of profiles: {len(profile_set)}")

The delta heuristic generates quartet profiles by examining distance patterns for 
each 4-taxon set.

Step 4: Infer Network
--------------------

Infer network from quartet profiles using SQuaRE algorithm:

.. code-block:: python

   from phylozoo.inference.squirrel import squirrel
   
   # Infer network from distance matrix (uses delta heuristic internally)
   network = squirrel(distances)
   
   print(f"Inferred network has {network.num_nodes} nodes")
   print(f"Level: {network.level()}")

The ``squirrel`` function combines delta heuristic and network reconstruction in one 
step. Alternatively, you can use the components separately:

.. code-block:: python

   from phylozoo.inference.squirrel import (
       delta_heuristic, quartet_joining, resolve_cycles
   )
   
   # Step-by-step inference
   profile_set = delta_heuristic(distances)
   tree = quartet_joining(profile_set)
   network = resolve_cycles(tree, profile_set)

Step 5: Validate and Analyze Network
-------------------------------------

Validate and analyze the inferred network:

.. code-block:: python

   from phylozoo.core.dnetwork import classifications, features
   
   # Check network properties
   print(f"Is tree: {network.is_tree()}")
   print(f"Level: {classifications.level(network)}")
   print(f"LSA node: {features.lsa_node(network)}")
   
   # Extract features
   blobs = features.blobs(network, trivial=False)
   print(f"Number of blobs: {len(blobs)}")

Step 6: Compare with Expected Network
--------------------------------------

If you have an expected network, compare quartet profiles:

.. code-block:: python

   from phylozoo.inference.squirrel import (
       sqprofileset_from_network, sqprofileset_similarity
   )
   from phylozoo import DirectedPhyNetwork
   
   # Load expected network
   expected_net = DirectedPhyNetwork.load("expected.enewick")
   
   # Extract quartet profiles from expected network
   expected_profiles = sqprofileset_from_network(expected_net)
   
   # Compare profiles
   similarity = sqprofileset_similarity(profile_set, expected_profiles)
   print(f"Profile similarity: {similarity}")

Step 7: Visualize Network
--------------------------

Visualize the inferred network:

.. code-block:: python

   from phylozoo.viz.dnetwork import plot_network
   
   # Plot network
   plot_network(network, show=True)
   
   # Save figure
   fig = plot_network(network, show=False)
   fig.savefig("inferred_network.png", dpi=300)

Step 8: Save Results
--------------------

Save the inferred network and intermediate results:

.. code-block:: python

   # Save network
   network.save("inferred_network.enewick")
   
   # Save distance matrix
   distances.save("distances.nexus")
   
   # Save quartet profiles (if format supported)
   # profile_set.save("profiles.txt")

Complete Example
----------------

Here's the complete inference workflow in one script:

.. code-block:: python

   from phylozoo import MSA, DistanceMatrix
   from phylozoo.core.sequence import hamming_distances
   from phylozoo.inference.squirrel import squirrel
   from phylozoo.viz.dnetwork import plot_network
   
   # Load sequences
   msa = MSA.load("alignment.fasta")
   
   # Compute distances
   distances = hamming_distances(msa)
   
   # Infer network
   network = squirrel(distances)
   
   # Analyze
   print(f"Inferred network level: {network.level()}")
   print(f"Number of hybrid nodes: {len(network.hybrid_nodes)}")
   
   # Visualize
   plot_network(network, show=True)
   
   # Save
   network.save("inferred_network.enewick")

Bootstrap Analysis
------------------

Assess confidence using bootstrapping:

.. code-block:: python

   from phylozoo.core.sequence import bootstrap
   from phylozoo.inference.squirrel import squirrel
   
   # Generate bootstrap replicates
   networks = []
   for bootstrapped_msa in bootstrap(msa, n_bootstrap=100, seed=42):
       bootstrapped_distances = hamming_distances(bootstrapped_msa)
       bootstrapped_network = squirrel(bootstrapped_distances)
       networks.append(bootstrapped_network)
   
   # Analyze bootstrap support
   # (e.g., count how often each edge appears)

.. note::
   The SQuaRE algorithm works best with accurate distance matrices. Consider 
   bootstrapping to assess confidence in inferred networks.

.. tip::
   Use the delta heuristic for quick quartet profile generation. For more accurate 
   profiles, consider using other methods or extracting profiles from gene trees.

.. seealso::
   For detailed information on inference functions, see :doc:`Inference <../manual/inference>`. 
   For network analysis workflows, see :doc:`Network Analysis Workflow <workflow_network_analysis>`.
