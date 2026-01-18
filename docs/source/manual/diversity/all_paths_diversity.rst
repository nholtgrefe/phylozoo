All-Paths Diversity
====================

The all-paths diversity measure computes phylogenetic diversity based on all paths from 
the root to a set of taxa in a phylogenetic network. It sums the branch lengths of all 
edges that lie on paths from the root to the selected taxa, providing a measure of the 
total evolutionary history captured by the taxon set.

Definition
----------

For a set of taxa :math:`S`, the all-paths diversity is the sum of branch lengths of all 
edges in the subnetwork induced by :math:`S` (i.e., all edges on paths from the root to 
any taxon in :math:`S`).

Mathematically, if :math:`E(S)` is the set of edges on paths from the root to taxa in 
:math:`S`, and :math:`\ell(e)` is the branch length of edge :math:`e`, then:

.. math::

   \text{All-Paths Diversity}(S) = \sum_{e \in E(S)} \ell(e)

This measure captures the total evolutionary history represented by the selected taxa, 
considering all possible paths through the network.

Basic Usage
-----------

Create an all-paths diversity measure and compute diversity:

.. code-block:: python

   from phylozoo.panda import AllPathsDiversity, diversity
   from phylozoo import DirectedPhyNetwork
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Create diversity measure
   measure = AllPathsDiversity()
   
   # Compute diversity for a set of taxa
   taxa_set = {"A", "B", "C", "D"}
   div_value = diversity(network, taxa_set, measure)
   
   print(f"Diversity: {div_value}")

Or use the convenience function:

.. code-block:: python

   from phylozoo.panda import all_paths
   
   # Create measure using convenience function
   measure = all_paths(network)
   
   # Compute diversity
   div_value = diversity(network, taxa_set, measure)

Marginal Diversity
------------------

Marginal diversity measures the change in total diversity when adding or removing a single 
taxon. For all-paths diversity, the marginal diversity of a taxon is the sum of branch 
lengths of edges that are unique to paths to that taxon (not shared with other selected taxa).

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

The greedy approach works for all networks and always produces a valid solution, though it 
may not be optimal.

Exact Maximization
------------------

For certain network structures, exact solutions can be found using dynamic programming 
based on scanwidth:

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

**Requirements for Exact Solution:**
* Network must be binary (or will be automatically converted to binary)
* Network must not have parallel edges
* Network must not have 2-blobs (blobs of size 2)

If these conditions are not met, the function raises ``PhyloZooNotImplementedError`` and 
you should use ``greedy_max_diversity`` instead.

Algorithm Details
-----------------

**Greedy Algorithm**
   The greedy algorithm iteratively selects the taxon with highest marginal diversity. 
   This is a simple heuristic that works for all networks but may not find the optimal solution.

**Exact Algorithm (Scanwidth-based DP)**
   The exact algorithm uses dynamic programming over a tree extension of the network. It 
   computes the scanwidth of the network and uses this to solve the maximum diversity 
   problem optimally. The algorithm has exponential time complexity in the scanwidth but 
   is fixed-parameter tractable.

**Binary Resolution**
   Non-binary networks are automatically converted to binary using binary resolution before 
   applying the exact algorithm. This ensures the algorithm works on all networks while 
   maintaining the diversity value.

Network Requirements
--------------------

**Greedy Maximization:**
* Works for all networks
* No restrictions on network structure
* Always produces a valid solution

**Exact Maximization:**
* Requires binary networks (automatically converted if needed)
* Requires networks without parallel edges
* Requires networks without 2-blobs
* Raises ``PhyloZooNotImplementedError`` if conditions are not met

Example: Complete Workflow
---------------------------

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.panda import (
       AllPathsDiversity, diversity, greedy_max_diversity, solve_max_diversity
   )
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Create measure
   measure = AllPathsDiversity()
   
   # Compute diversity for specific taxa
   taxa = {"A", "B", "C"}
   div = diversity(network, taxa, measure)
   print(f"Diversity of {taxa}: {div}")
   
   # Try exact solution first
   k = 5
   try:
       opt_div, opt_taxa = solve_max_diversity(network, k, measure)
       print(f"Optimal {k} taxa (exact): {opt_taxa}")
       print(f"Optimal diversity: {opt_div}")
   except PhyloZooNotImplementedError:
       # Fall back to greedy if exact not available
       opt_div, opt_taxa = greedy_max_diversity(network, k, measure)
       print(f"Optimal {k} taxa (greedy): {opt_taxa}")
       print(f"Diversity: {opt_div}")

Comparison with Other Measures
------------------------------

All-paths diversity differs from other diversity measures in that it considers **all** paths 
from the root to the selected taxa, not just a single path per taxon. This makes it 
particularly suitable for networks with reticulations, where multiple paths may exist.

**Advantages:**
* Captures full evolutionary history including reticulations
* Works well with hybrid nodes and cycles
* Can be computed exactly for certain network structures

**Considerations:**
* Exact solution requires specific network structures
* May overcount shared evolutionary history in some cases
* Computational complexity increases with network complexity

Best Practices
--------------

1. **Use greedy for quick results**: When you need a solution quickly or your network doesn't 
   meet the requirements for exact solution, use ``greedy_max_diversity``.

2. **Try exact when possible**: If your network meets the requirements, use ``solve_max_diversity`` 
   for optimal solutions.

3. **Handle exceptions**: Always handle ``PhyloZooNotImplementedError`` when using exact 
   maximization, and fall back to greedy if needed.

4. **Validate networks**: Ensure your network has branch lengths on all edges, as the measure 
   depends on branch length values.

5. **Consider network structure**: The exact algorithm works best on tree-like networks or 
   networks with low scanwidth. For highly complex networks, greedy may be more practical.

.. seealso::
   For more information on:
   * Diversity framework: :doc:`Diversity Module <../diversity/index>`
   * Network structures: :doc:`Networks (Advanced) <../core/networks/advanced>`
   * Scanwidth: See the implementation in ``phylozoo.panda.utils.scanwidth``
