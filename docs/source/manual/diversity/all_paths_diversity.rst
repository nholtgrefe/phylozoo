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
   
   measure = all_paths(network)
   div_value = diversity(network, taxa_set, measure)

Maximization
-----------

**Greedy Maximization**

Greedy maximization iteratively selects taxa to maximize marginal diversity. It works 
for all networks and always produces a valid solution, though it may not be optimal:

.. code-block:: python

   from phylozoo.panda import greedy_max_diversity
   
   k = 5
   diversity_value, selected_taxa = greedy_max_diversity(network, k, measure)
   print(f"Selected taxa: {selected_taxa}")
   print(f"Diversity value: {diversity_value}")

**Exact Maximization**

For certain network structures, exact solutions can be found using dynamic programming 
based on scanwidth :cite:`Holtgrefe2024Exact`:

.. code-block:: python

   from phylozoo.panda import solve_max_diversity
   
   k = 5
   try:
       diversity_value, selected_taxa = solve_max_diversity(network, k, measure)
       print(f"Optimal taxa: {selected_taxa}")
       print(f"Optimal diversity: {diversity_value}")
   except PhyloZooNotImplementedError:
       # Fall back to greedy if exact not available
       diversity_value, selected_taxa = greedy_max_diversity(network, k, measure)

**Requirements for Exact Solution:**
* Network must be binary (automatically converted if needed)
* Network must not have parallel edges
* Network must not have 2-blobs (blobs of size 2)

If these conditions are not met, the function raises ``PhyloZooNotImplementedError``.

Algorithm Details
-----------------

The exact algorithm uses dynamic programming over a tree extension of the network, 
computing the scanwidth :cite:`Holtgrefe2024Exact` to solve the maximum diversity problem 
optimally. The algorithm has exponential time complexity in the scanwidth but is 
fixed-parameter tractable. Non-binary networks are automatically converted to binary using 
binary resolution before applying the exact algorithm.

.. note::
   The scanwidth computation will be moved to a separate package in a future version of 
   PhyloZoo, which will become a new dependency. This will allow for better modularity and 
   independent development of scanwidth algorithms.

Example: Complete Workflow
---------------------------

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.panda import (
       AllPathsDiversity, diversity, greedy_max_diversity, solve_max_diversity
   )
   
   network = DirectedPhyNetwork.load("network.enewick")
   measure = AllPathsDiversity()
   
   # Compute diversity for specific taxa
   taxa = {"A", "B", "C"}
   div = diversity(network, taxa, measure)
   print(f"Diversity of {taxa}: {div}")
   
   # Try exact solution first, fall back to greedy
   k = 5
   try:
       opt_div, opt_taxa = solve_max_diversity(network, k, measure)
       print(f"Optimal {k} taxa (exact): {opt_taxa}")
   except PhyloZooNotImplementedError:
       opt_div, opt_taxa = greedy_max_diversity(network, k, measure)
       print(f"Optimal {k} taxa (greedy): {opt_taxa}")

Best Practices
--------------

1. **Use greedy for quick results**: When you need a solution quickly or your network doesn't 
   meet the requirements for exact solution, use ``greedy_max_diversity``.

2. **Try exact when possible**: If your network meets the requirements, use ``solve_max_diversity`` 
   for optimal solutions. Always handle ``PhyloZooNotImplementedError`` and fall back to greedy.

3. **Ensure branch lengths**: The measure depends on branch length values. Ensure your network 
   has branch lengths on all edges.

4. **Consider network structure**: The exact algorithm works best on tree-like networks or 
   networks with low scanwidth. For highly complex networks, greedy may be more practical.

.. seealso::
   For more information on:
   * Diversity framework: :doc:`Diversity Module <../diversity/index>`
   * Network structures: :doc:`Networks (Advanced) <../core/networks/advanced>`
