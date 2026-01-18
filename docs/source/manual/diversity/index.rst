Diversity Module (Panda)
=========================

The diversity module (`phylozoo.panda`) provides a flexible framework for calculating 
phylogenetic diversity and finding optimal sets of taxa that maximize diversity 
:cite:`PhyloZoo2024`.

For detailed information on the all-paths diversity measure, see :doc:`All-Paths Diversity <all_paths_diversity>`.

Diversity Measures
------------------

PhyloZoo uses a protocol-based architecture for diversity measures. The base protocol 
is ``DiversityMeasure``, which defines the interface that all diversity measures must 
implement.

Currently implemented measures:

* **All-Paths Diversity**: Measures diversity based on all paths between taxa in the 
  network. Considers all paths from root to each taxon, weighted by edge probabilities. 
  See :doc:`All-Paths Diversity <all_paths_diversity>` for details.

Creating a Diversity Measure
-----------------------------

.. code-block:: python

   from phylozoo.panda import AllPathsDiversity
   from phylozoo import DirectedPhyNetwork
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Create diversity measure
   measure = AllPathsDiversity(network)

Computing Diversity
-------------------

Compute diversity for a set of taxa:

.. code-block:: python

   from phylozoo.panda import diversity
   
   # Compute diversity for a set of taxa
   taxa_set = {"A", "B", "C", "D"}
   div_value = diversity(network, taxa_set, measure)
   
   print(f"Diversity: {div_value}")

Marginal Diversity
------------------

Marginal diversity measures the change in total diversity when adding or removing a 
single taxon:

.. code-block:: python

   from phylozoo.panda import marginal_diversities
   
   # Compute marginal diversity for all taxa
   current_set = {"A", "B"}
   marginals = marginal_diversities(network, current_set, measure)
   
   # marginals is a dict mapping taxon -> marginal diversity value
   for taxon, marginal in marginals.items():
       print(f"{taxon}: {marginal}")

Greedy Maximization
-------------------

Greedy maximization iteratively selects taxa to maximize marginal diversity:

.. code-block:: python

   from phylozoo.panda import greedy_max_diversity
   
   # Find k taxa that maximize diversity (greedy approach)
   k = 5
   diversity_value, selected_taxa = greedy_max_diversity(network, k, measure)
   
   print(f"Selected taxa: {selected_taxa}")
   print(f"Diversity value: {diversity_value}")

The greedy approach works for all networks and diversity measures.

Exact Maximization
------------------

For certain network structures, exact solutions can be found using dynamic programming:

.. code-block:: python

   from phylozoo.panda import solve_max_diversity
   
   # Find k taxa that maximize diversity (exact solution)
   k = 5
   try:
       diversity_value, selected_taxa = solve_max_diversity(network, k, measure)
       print(f"Optimal taxa: {selected_taxa}")
       print(f"Optimal diversity: {diversity_value}")
   except PhyloZooNotImplementedError:
       print("Exact solution not available, use greedy_max_diversity")

The exact solution method is only available for certain diversity measures and network 
structures (e.g., networks without 2-blobs or parallel edges).

Complete Example
----------------

Here's a complete example:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.panda import (
       AllPathsDiversity, diversity, greedy_max_diversity
   )
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Create measure
   measure = AllPathsDiversity(network)
   
   # Compute diversity for specific taxa
   taxa = {"A", "B", "C"}
   div = diversity(network, taxa, measure)
   
   # Find optimal k taxa
   k = 5
   opt_div, opt_taxa = greedy_max_diversity(network, k, measure)
   
   print(f"Optimal {k} taxa: {opt_taxa}")
   print(f"Diversity: {opt_div}")

Available Functions
-------------------

**Diversity Computation:**

* **diversity(network, taxa, measure, **kwargs)** - Compute diversity of a set of taxa. 
  Validates that all taxa are in the network. Returns float diversity value. See 
  :func:`phylozoo.panda.diversity`.

* **marginal_diversities(network, saved_taxa, measure, **kwargs)** - Compute marginal 
  diversity contributions for all taxa. For taxa in saved_taxa, returns decrease when 
  removing. For taxa not in saved_taxa, returns increase when adding. Returns dict 
  mapping taxon -> marginal diversity. See :func:`phylozoo.panda.marginal_diversities`.

**Optimization:**

* **greedy_max_diversity(network, k, measure, **kwargs)** - Greedily find k taxa that 
  maximize diversity. Iteratively selects taxa with highest marginal diversity. Works 
  for all networks and measures. Returns tuple (diversity_value, set_of_taxa). See 
  :func:`phylozoo.panda.greedy_max_diversity`.

* **solve_max_diversity(network, k, measure, **kwargs)** - Solve maximum diversity 
  problem exactly using dynamic programming. Only available for certain measures and 
  network structures. Raises PhyloZooNotImplementedError if not available. Returns 
  tuple (diversity_value, set_of_taxa). See :func:`phylozoo.panda.solve_max_diversity`.

**Diversity Measures:**

* **AllPathsDiversity(network)** - All-paths diversity measure. Measures diversity 
  based on all paths from root to taxa, weighted by edge probabilities. Requires 
  networks without parallel edges and 2-blobs for exact solutions. See 
  :class:`phylozoo.panda.AllPathsDiversity`.

* **all_paths(network)** - Convenience function to create AllPathsDiversity measure. 
  Returns AllPathsDiversity instance. See :func:`phylozoo.panda.all_paths`.

**Protocol:**

* **DiversityMeasure** - Protocol defining interface for diversity measures. All 
  measures must implement ``compute_diversity`` method. See 
  :class:`phylozoo.panda.DiversityMeasure`.

Network Requirements
--------------------

Different diversity measures have different requirements:

* **All-Paths Diversity**: 
  - Greedy maximization: Works for all networks
  - Exact maximization: Requires networks without parallel edges and without 2-blobs

.. note::
   The exact solution method uses scanwidth-based dynamic programming, which requires 
   specific network structures. If exact solution is not available, use greedy 
   maximization instead.

.. tip::
   Use greedy maximization for quick results on any network. Use exact maximization 
   when you need optimal solutions and your network structure supports it.

.. warning::
   The ``solve_max_diversity`` function will raise ``PhyloZooNotImplementedError`` if 
   exact solution is not available. Always handle this exception or use greedy 
   maximization as a fallback.

.. toctree::
   :maxdepth: 1
   :hidden:
   
   all_paths_diversity

.. seealso::
   For detailed measure information, see :doc:`All-Paths Diversity <all_paths_diversity>`.
   For network analysis workflows, see :doc:`Network Analysis Workflow <../../tutorials/workflow_network_analysis>`. 
   For I/O operations, see :doc:`I/O <../io>`.
