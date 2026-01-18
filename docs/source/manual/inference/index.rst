Inference Module
================

The inference module (`phylozoo.inference`) provides algorithms for inferring 
phylogenetic networks from data :cite:`PhyloZoo2024`. The main algorithm is SQuaRE 
(SQuirrel quartet-based REconstruction), which infers semi-directed networks from 
quartet profiles.

For detailed information on the SQuaRE algorithm, see :doc:`SQuaRE Algorithm <squirrel>`.

SQuaRE Algorithm Overview
--------------------------

The SQuaRE algorithm reconstructs phylogenetic networks from quartet profiles using 
a multi-step process:

1. **Quartet Profile Generation**: Generate quartet profiles from input data (e.g., 
   distance matrices)
2. **Quartet Joining**: Build an initial tree structure using quartet joining
3. **Cycle Resolution**: Identify and resolve cycles to add reticulations
4. **Network Assembly**: Combine components into a complete network

Basic Usage
-----------

The simplest way to infer a network is using the ``squirrel`` function:

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

Step-by-Step Inference
-----------------------

For more control, use the components separately:

.. code-block:: python

   from phylozoo.inference.squirrel import (
       delta_heuristic, adapted_quartet_joining, resolve_cycles
   )
   
   # Step 1: Generate quartet profiles
   profile_set = delta_heuristic(distance_matrix)
   
   # Step 2: Build initial tree
   tree = adapted_quartet_joining(profile_set)
   
   # Step 3: Resolve cycles to add reticulations
   network = resolve_cycles(tree, profile_set)

Complete Workflow Example
-------------------------

See :doc:`Inference Workflow <../../tutorials/workflow_inference>` for a complete example showing 
the full inference pipeline from sequences to networks.

Available Functions
-------------------

**Main Inference Function:**

* **squirrel(profileset, outgroup=None, **kwargs)** - Main SQuaRE algorithm. Reconstructs 
  phylogenetic network from squirrel quartet profile set. Can take DistanceMatrix directly 
  (uses delta heuristic internally) or SqQuartetProfileSet. Returns SemiDirectedPhyNetwork 
  or DirectedPhyNetwork (if outgroup specified). See :func:`phylozoo.inference.squirrel.squirrel`.

**Quartet Profile Generation:**

* **delta_heuristic(distance_matrix, lam=0.3, weight=True)** - Generate quartet profiles 
  from distance matrix using delta heuristic. Uses distance patterns to determine if quartets 
  should be resolved as splits or cycles. Returns SqQuartetProfileSet. See 
  :func:`phylozoo.inference.squirrel.delta_heuristic`.

* **sqprofileset_from_network(network)** - Extract quartet profiles from existing network. 
  Useful for comparing networks or generating expected quartet distributions. Returns 
  SqQuartetProfileSet. See :func:`phylozoo.inference.squirrel.sqprofileset_from_network`.

**Quartet Joining:**

* **quartet_joining(profileset)** - Standard quartet joining algorithm. Builds tree from 
  quartet profiles. Requires fully resolved quartets. Returns SemiDirectedPhyNetwork. See 
  :func:`phylozoo.inference.squirrel.quartet_joining`.

* **adapted_quartet_joining(profileset, starting_tree=None)** - Adapted quartet joining 
  that handles unresolved quartets. Can use optional starting tree. Returns 
  SemiDirectedPhyNetwork. See :func:`phylozoo.inference.squirrel.adapted_quartet_joining`.

**Tree Building:**

* **tstar_tree(profileset)** - Build T* tree from quartet profiles. T* is a tree that 
  maximizes quartet support. Returns SemiDirectedPhyNetwork. See 
  :func:`phylozoo.inference.squirrel.tstar_tree`.

* **bstar(profileset)** - Build B* split system from quartet profiles. Returns SplitSystem. 
  See :func:`phylozoo.inference.squirrel.bstar`.

**Cycle Resolution:**

* **resolve_cycles(network, profileset, rho=(0.5, 1.0, 0.5, 1.0), tsp_threshold=13)** - 
  Resolve cycles in network to add reticulations. Uses quartet profiles to determine cycle 
  structure. Returns SemiDirectedPhyNetwork. See 
  :func:`phylozoo.inference.squirrel.resolve_cycles`.

**Tree Unresolution:**

* **unresolve_tree(tree, profileset)** - Iteratively contract least-supported splits to 
  generate contraction sequence. Yields SemiDirectedPhyNetwork objects. See 
  :func:`phylozoo.inference.squirrel.unresolve_tree`.

* **split_support(tree, profileset)** - Compute support for each split in tree based on 
  quartet profiles. Returns dictionary mapping splits to support values. See 
  :func:`phylozoo.inference.squirrel.split_support`.

**Comparison:**

* **sqprofileset_similarity(profileset1, profileset2)** - Compute similarity between two 
  quartet profile sets. Returns float similarity score. See 
  :func:`phylozoo.inference.squirrel.sqprofileset_similarity`.

**Utility Functions:**

* **root_at_outgroup(network, outgroup)** - Root network at outgroup taxon. Converts 
  SemiDirectedPhyNetwork to DirectedPhyNetwork. Returns DirectedPhyNetwork. See 
  :func:`phylozoo.inference.utils.root_at_outgroup`.

**Classes:**

* **SqQuartetProfile** - Squirrel quartet profile. Contains 1 or 2 resolved quartets, 
  optionally with reticulation leaf. See :class:`phylozoo.inference.squirrel.SqQuartetProfile`.

* **SqQuartetProfileSet** - Collection of squirrel quartet profiles. Must be dense (profile 
  for every 4-taxon set). See :class:`phylozoo.inference.squirrel.SqQuartetProfileSet`.

.. note::
   The SQuaRE algorithm works best with accurate distance matrices. Consider bootstrapping 
   to assess confidence in inferred networks.

.. tip::
   Use ``squirrel`` with a DistanceMatrix for the simplest workflow. For more control, use 
   the components separately.

.. warning::
   The delta heuristic requires a dense distance matrix (distances for all pairs). Missing 
   values will cause errors.

.. toctree::
   :maxdepth: 1
   :hidden:
   
   squirrel

.. seealso::
   For detailed algorithm information, see :doc:`SQuaRE Algorithm <squirrel>`.
   For a complete inference workflow example, see :doc:`Inference Workflow <../../tutorials/workflow_inference>`. 
   For network analysis workflows, see :doc:`Network Analysis Workflow <../../tutorials/workflow_network_analysis>`.
