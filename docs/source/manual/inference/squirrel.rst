SQuaRE Algorithm (Squirrel)
===========================

The SQuaRE (SQuirrel quartet-based REconstruction) algorithm is the main network inference 
method in PhyloZoo. It reconstructs semi-directed phylogenetic networks from quartet profiles 
using a multi-step process that combines quartet joining with cycle resolution.

Algorithm Overview
------------------

The SQuaRE algorithm works in four main steps:

1. **T* Tree Construction**: Builds a tree that maximizes quartet support from the profile set
2. **Quartet Joining**: Applies adapted quartet joining to build an initial tree structure
3. **Tree Unresolution**: Iteratively contracts least-supported splits to generate a sequence 
   of contracted trees
4. **Cycle Resolution**: For each contracted tree, resolves cycles to add reticulations and 
   computes similarity scores, returning the network with highest similarity

The algorithm explores multiple network topologies by systematically contracting splits and 
evaluating which network best matches the input quartet profiles.

Basic Usage
-----------

The simplest way to use SQuaRE is with the ``squirrel`` function:

.. code-block:: python

   from phylozoo.inference.squirrel import squirrel
   from phylozoo import DistanceMatrix
   
   # Infer network from distance matrix (uses delta heuristic internally)
   distance_matrix = DistanceMatrix.load("distances.nexus")
   network = squirrel(distance_matrix)
   
   # The result is a SemiDirectedPhyNetwork
   print(f"Inferred network level: {network.level()}")

You can also specify an outgroup to get a rooted network:

.. code-block:: python

   # Infer rooted network
   rooted_network = squirrel(distance_matrix, outgroup="outgroup_taxon")
   # Returns DirectedPhyNetwork

Step-by-Step Process
--------------------

For more control, you can use the components separately:

.. code-block:: python

   from phylozoo.inference.squirrel import (
       delta_heuristic, adapted_quartet_joining, unresolve_tree, resolve_cycles
   )
   
   # Step 1: Generate quartet profiles from distance matrix
   profile_set = delta_heuristic(distance_matrix)
   
   # Step 2: Build initial tree using quartet joining
   tree = adapted_quartet_joining(profile_set)
   
   # Step 3: Iterate through contracted trees
   for contracted_tree in unresolve_tree(tree, profile_set):
       # Step 4: Resolve cycles to add reticulations
       network = resolve_cycles(contracted_tree, profile_set)
       # Evaluate network...

The ``squirrel`` function automates this process and selects the best network.

Parallelization
---------------

SQuaRE supports parallelization to speed up the evaluation of multiple contracted trees:

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

**T* Tree**
   The T* tree is a tree that maximizes quartet support. It serves as a starting point for 
   quartet joining. Computed using :func:`phylozoo.inference.squirrel.tstar_tree`.

**Quartet Joining**
   Adapted quartet joining builds a tree structure from quartet profiles, handling both 
   resolved and unresolved quartets. Uses :func:`phylozoo.inference.squirrel.adapted_quartet_joining`.

**Tree Unresolution**
   Iteratively contracts least-supported splits to generate a sequence of contracted trees. 
   Each contraction removes one split, creating progressively simpler topologies. Uses 
   :func:`phylozoo.inference.squirrel.unresolve_tree`.

**Cycle Resolution**
   For each contracted tree, identifies cycles in the quartet profiles and resolves them by 
   adding reticulations. Computes similarity scores to evaluate how well the network matches 
   the profiles. Uses :func:`phylozoo.inference.squirrel.resolve_cycles`.

**Similarity Scoring**
   The algorithm evaluates each candidate network by computing its similarity to the input 
   quartet profiles. The network with the highest similarity score is returned. Uses 
   :func:`phylozoo.inference.squirrel.sqprofileset_similarity`.

Parameters
----------

The ``squirrel`` function accepts several parameters:

* **profileset**: ``SqQuartetProfileSet`` - The quartet profile set (must be dense)
* **outgroup**: ``str | None`` - Optional outgroup taxon for rooting
* **parallel**: ``ParallelConfig | None`` - Parallelization configuration
* **rho**: ``tuple[float, float, float, float]`` - Rho vector for distance computation in 
  cycle resolution (default: ``(0.5, 1.0, 0.5, 1.0)``)
* **tsp_threshold**: ``int | None`` - Maximum partition size for optimal TSP solving 
  (default: ``13``)

Example: Complete Workflow
---------------------------

.. code-block:: python

   from phylozoo import DistanceMatrix
   from phylozoo.inference.squirrel import squirrel
   from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
   
   # Load distance matrix
   dm = DistanceMatrix.load("distances.nexus")
   
   # Infer network with parallelization
   network = squirrel(
       dm,  # Can pass DistanceMatrix directly
       outgroup="outgroup_taxon",  # Optional: root the network
       parallel=ParallelConfig(
           backend=ParallelBackend.MULTIPROCESSING,
           n_jobs=4
       )
   )
   
   # Save inferred network
   network.save("inferred_network.newick")
   print(f"Network level: {network.level()}")
   print(f"Number of hybrid nodes: {len(network.hybrid_nodes())}")

Limitations and Considerations
------------------------------

* **Dense Profile Sets**: The algorithm requires a dense quartet profile set (a profile for 
  every 4-taxon set). Use ``delta_heuristic`` to generate profiles from distance matrices.

* **Computational Complexity**: The algorithm explores multiple network topologies, which 
  can be computationally intensive for large datasets. Parallelization helps but may not 
  scale linearly.

* **Quality of Input**: The quality of inferred networks depends on the quality of input 
  quartet profiles. Accurate distance matrices lead to better quartet profiles and thus 
  better networks.

* **Network Level**: The algorithm can infer networks of various levels, but the level 
  depends on the quartet profile structure and the contraction sequence.

Best Practices
--------------

1. **Use accurate distance matrices**: The quality of inferred networks depends on accurate 
   input data. Consider bootstrapping to assess confidence.

2. **Use parallelization for large datasets**: When processing many taxa or complex networks, 
   parallelization can significantly speed up inference.

3. **Experiment with parameters**: The ``rho`` vector and ``tsp_threshold`` can affect 
   results. Experiment with different values if needed.

4. **Validate results**: Always validate inferred networks and consider comparing multiple 
   runs or using different starting trees.

.. seealso::
   For more information on:
   * Quartet profiles: :doc:`Quartets <../core/quartets>`
   * Distance matrices: :doc:`Distance Matrices <../core/distance>`
   * Parallelization: :doc:`Parallelization <../advanced/parallelization>`
   * Complete inference workflow: :doc:`Inference Workflow <../../tutorials/workflow_inference>`
