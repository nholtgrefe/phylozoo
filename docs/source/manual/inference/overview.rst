Overview
========

The inference module (`phylozoo.inference`) provides algorithms for inferring 
phylogenetic networks from biological data. This module contains methods for 
reconstructing evolutionary relationships from distance matrices, sequence alignments, 
and other phylogenetic data sources.

Available Algorithms
--------------------

**Squirrel (Semi-directed Quarnet-based Inference to Reconstruct Level-1 Networks)**
   The main network inference algorithm in PhyloZoo. Reconstructs semi-directed 
   phylogenetic level-1 networks from quartet profiles using a multi-step process 
   :cite:`Holtgrefe2025Squirrel`. Quartet profiles are generated from distance matrices 
   using the delta heuristic. Distances can be computed from sequence alignments using 
   Hamming distance.

   * **Input**: Quartet profiles (typically generated from distance matrices)
   * **Output**: Semi-directed or directed phylogenetic level-1 networks
   * **Key Features**: Handles reticulate evolution, parallelizable, explores multiple 
     network topologies

   For detailed information, see :doc:`Squirrel Algorithm <squirrel>`.

General Workflow
----------------

Most inference workflows in PhyloZoo follow this pattern:

1. **Data Preparation**: Load biological input data (sequence alignments, distance matrices)
2. **Feature Extraction**: Convert data into quartet profiles (e.g., using delta heuristic)
3. **Network Inference**: Apply inference algorithms (e.g., Squirrel)
4. **Post-processing**: Root networks, validate results, analyze properties

.. code-block:: python

   from phylozoo import DistanceMatrix
   from phylozoo.inference.squirrel import delta_heuristic, squirrel
   
   # 1. Load data
   distances = DistanceMatrix.load("input_distances.nexus")
   
   # 2. Extract features
   quartet_profiles = delta_heuristic(distances)
   
   # 3. Infer network
   network = squirrel(quartet_profiles)
   
   # 4. Post-process
   rooted_network = squirrel(quartet_profiles, outgroup="outgroup_taxon")

.. seealso::
   For complete inference workflow examples, see :doc:`Inference Workflow <../../tutorials/workflow_inference>`.
   For network analysis after inference, see :doc:`Network Analysis Workflow <../../tutorials/workflow_network_analysis>`.
