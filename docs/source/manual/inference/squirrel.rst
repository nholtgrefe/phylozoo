Squirrel Algorithm
===========================

The Squirrel (Semi-directed Quarnet-based Inference to Reconstruct Level-1 Networks) algorithm is the main network inference 
method in PhyloZoo. It reconstructs semi-directed phylogenetic networks from quartet profiles 
using a multi-step process :cite:`Holtgrefe2025Squirrel`. The quartet profiles are generated from distance matrices using the delta heuristic.

Algorithm Overview
------------------

The Squirrel algorithm works in six main steps:

1. **Quartet Profile Generation**: Generate quartet profiles from distance matrices using the delta heuristic
2. **T* Tree Construction**: Builds a tree that maximizes quartet support from the profile set
3. **Quartet Joining**: Applies adapted quartet joining to build an initial tree structure
4. **Tree Unresolution**: Iteratively contracts least-supported splits to generate a sequence 
   of contracted trees
5. **Cycle Resolution**: For each contracted tree, resolves cycles to add reticulations and 
   computes similarity scores
6. **Network Selection**: Returns the network with the highest similarity score

The algorithm explores multiple network topologies by systematically contracting splits and 
evaluating which network best matches the input quartet profiles.

Basic Usage
-----------

The simplest way to use Squirrel is with the ``squirrel`` function:

.. code-block:: python

   from phylozoo.inference.squirrel import squirrel, delta_heuristic
   from phylozoo import DistanceMatrix
   
   # Load distance matrix and generate quartet profiles
   distance_matrix = DistanceMatrix.load("distances.nexus")
   profile_set = delta_heuristic(distance_matrix)
   
   # Infer network from quartet profile set
   network = squirrel(profile_set)
   
   # The result is a SemiDirectedPhyNetwork
   print(f"Inferred network level: {network.level()}")

You can also specify an outgroup to get a rooted network:

.. code-block:: python

   # Infer rooted network
   rooted_network = squirrel(profile_set, outgroup="outgroup_taxon")
   # Returns DirectedPhyNetwork

.. tip::
   **Computing Distances from Sequence Alignments**
   
   If you start with sequence alignments instead of distance matrices, you can compute 
   distances using the ``hamming_distances`` function:
   
   .. code-block:: python
   
      from phylozoo import MSA
      from phylozoo.core.sequence import hamming_distances
      from phylozoo.inference.squirrel import delta_heuristic, squirrel
      
      # Load sequence alignment
      msa = MSA.load("alignment.fasta")
      
      # Compute distance matrix from sequences
      distance_matrix = hamming_distances(msa)
      
      # Generate quartet profiles and infer network
      profile_set = delta_heuristic(distance_matrix)
      network = squirrel(profile_set)
   
   The ``hamming_distances`` function computes normalized Hamming distances, excluding 
   positions with gaps or unknown characters. See :doc:`Sequences <../core/sequences>` 
   for more details on sequence operations.

Step-by-Step Process
--------------------

For more control, you can use the components separately:

.. code-block:: python

   from phylozoo.inference.squirrel import (
       delta_heuristic, tstar_tree, adapted_quartet_joining, 
       unresolve_tree, resolve_cycles
   )
   
   # Step 1: Generate quartet profiles from distance matrix
   profile_set = delta_heuristic(distance_matrix)
   
   # Step 2: Build T* tree (maximizes quartet support)
   tstar = tstar_tree(profile_set)
   
   # Step 3: Build initial tree using quartet joining
   tree = adapted_quartet_joining(profile_set, starting_tree=tstar)
   
   # Step 4: Iterate through contracted trees
   for contracted_tree in unresolve_tree(tree, profile_set):
       # Step 5: Resolve cycles to add reticulations
       network = resolve_cycles(profile_set, contracted_tree)
       # Evaluate network...

The ``squirrel`` function automates this process and selects the best network.

Parallelization
---------------

Squirrel supports parallelization to speed up the evaluation of multiple contracted trees:

.. code-block:: python

   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   from phylozoo.inference.squirrel import squirrel
   
   # Use multiprocessing with 4 cores
   network = squirrel(
       profileset,
       parallel=ParallelConfig(
           backend=ParallelBackend.MULTIPROCESSING,
           n_jobs=4
       )
   )
   
   # Use all available cores
   network = squirrel(
       profileset,
       parallel=ParallelConfig(
           backend=ParallelBackend.MULTIPROCESSING,
           n_jobs=None  # auto-detect
       )
   )

Parallelization is most beneficial when processing many contracted trees, as each tree 
can be evaluated independently.

Algorithm Components
--------------------

The algorithm consists of several components that work together:

* **Quartet Profile Generation**: Converts distance matrices to quartet profiles using the 
  delta heuristic (:func:`phylozoo.inference.squirrel.delta_heuristic`)
* **T* Tree Construction**: Builds a tree that maximizes quartet support 
  (:func:`phylozoo.inference.squirrel.tstar_tree`)
* **Quartet Joining**: Builds an initial tree structure from quartet profiles 
  (:func:`phylozoo.inference.squirrel.adapted_quartet_joining`)
* **Tree Unresolution**: Iteratively contracts least-supported splits to generate contracted 
  trees (:func:`phylozoo.inference.squirrel.unresolve_tree`)
* **Cycle Resolution**: Resolves cycles to add reticulations and computes similarity scores 
  (:func:`phylozoo.inference.squirrel.resolve_cycles`)
* **Network Selection**: Returns the network with the highest similarity score

Parameters
----------

The ``squirrel`` function accepts:

* **profileset**: ``SqQuartetProfileSet`` - The quartet profile set (must be dense)
* **outgroup**: ``str | None`` - Optional outgroup taxon for rooting
* **parallel**: ``ParallelConfig | None`` - Parallelization configuration
* **rho**: ``tuple[float, float, float, float]`` - Rho vector for distance computation in 
  cycle resolution (default: ``(0.5, 1.0, 0.5, 1.0)``)
* **tsp_threshold**: ``int | None`` - Maximum partition size for optimal TSP solving 
  (default: ``13``)

For detailed API documentation of all functions and classes, see the 
:mod:`phylozoo.inference.squirrel` module reference.

Example: Complete Workflow
---------------------------

.. code-block:: python

   from phylozoo import DistanceMatrix
   from phylozoo.inference.squirrel import squirrel, delta_heuristic
   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   
   # Load distance matrix and generate quartet profiles
   dm = DistanceMatrix.load("distances.nexus")
   profile_set = delta_heuristic(dm)
   
   # Infer network with parallelization
   network = squirrel(
       profile_set,
       outgroup="outgroup_taxon",  # Optional: root the network
       parallel=ParallelConfig(
           backend=ParallelBackend.MULTIPROCESSING,
           n_jobs=4
       )
   )
   
   # Save and analyze
   network.save("inferred_network.newick")
   print(f"Network level: {network.level()}")
   print(f"Number of hybrid nodes: {len(network.hybrid_nodes())}")

Limitations and Best Practices
-------------------------------

* **Dense Profile Sets**: The algorithm requires a dense quartet profile set (a profile for 
  every 4-taxon set). Use ``delta_heuristic`` to generate profiles from distance matrices.

* **Computational Complexity**: The algorithm explores multiple network topologies, which 
  can be computationally intensive for large datasets. Use parallelization for large datasets.

* **Quality of Input**: The quality of inferred networks depends on accurate input data. 
  Consider bootstrapping to assess confidence.

* **Parameter Tuning**: The ``rho`` vector and ``tsp_threshold`` parameters can affect 
  results. Experiment with different values if needed.

* **Validation**: Always validate inferred networks and consider comparing multiple runs.

.. warning::
   The delta heuristic requires a dense distance matrix (distances for all pairs). Missing 
   values will cause errors.

.. seealso::
   For more information on:
   * Quartet profiles: :doc:`Quartets <../core/quartets>`
   * Distance matrices: :doc:`Distance Matrices <../core/distance>`
   * Parallelization: :doc:`Parallelization <../advanced/parallelization>`
   * Complete inference workflow: :doc:`Inference Workflow <../../tutorials/workflow_inference>`
