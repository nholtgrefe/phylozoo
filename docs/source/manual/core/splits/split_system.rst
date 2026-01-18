SplitSystem
===========

A ``SplitSystem`` is a collection of splits where each split covers the complete set 
of elements. It validates that all splits cover the same element set and supports I/O 
operations.

Creating Split Systems
----------------------

.. code-block:: python

   from phylozoo import Split, SplitSystem
   
   splits = [
       Split({"A", "B"}, {"C", "D"}),
       Split({"A"}, {"B", "C", "D"}),
       Split({"A", "B", "C"}, {"D"})
   ]
   
   split_system = SplitSystem(splits)
   
   # Access properties
   elements = split_system.elements
   num_splits = len(split_system)

Weighted Split Systems
----------------------

A ``WeightedSplitSystem`` associates weights with splits:

.. code-block:: python

   from phylozoo.core.split import WeightedSplitSystem
   
   weighted_system = WeightedSplitSystem(
       splits=splits,
       weights=[0.9, 0.8, 0.7],
       elements={"A", "B", "C", "D"}
   )
   
   # Access weights
   weight = weighted_system.get_weight(splits[0])
   
   # Convert from regular split system
   from phylozoo.core.split import to_weightedsplitsystem
   weighted = to_weightedsplitsystem(split_system, default_weight=1.0)

Split System Algorithms
-----------------------

Convert split systems to trees and extract other structures:

.. code-block:: python

   from phylozoo.core.split import algorithms
   
   # Convert split system to tree (if compatible)
   tree = algorithms.splitsystem_to_tree(split_system)
   
   # Extract distance matrix from split system
   distances = algorithms.distances_from_splitsystem(split_system)
   
   # Extract quartets from split system
   quartets = algorithms.quartets_from_splitsystem(split_system)

Complete Example
----------------

.. code-block:: python

   from phylozoo import Split, SplitSystem
   from phylozoo.core.split import algorithms, classifications
   
   # Create splits representing phylogenetic relationships
   splits = [
       Split({"A", "B"}, {"C", "D", "E"}),  # A,B vs C,D,E
       Split({"A"}, {"B", "C", "D", "E"}),  # A vs others
       Split({"A", "B", "C"}, {"D", "E"}),  # A,B,C vs D,E
   ]
   
   # Create split system
   system = SplitSystem(splits)
   
   # Check if splits are tree-compatible
   compatible = classifications.is_tree_compatible(system)
   
   if compatible:
       # Convert to tree
       tree = algorithms.splitsystem_to_tree(system)
       
       # Extract distance matrix
       distances = algorithms.distances_from_splitsystem(system)
       
       # Extract quartets
       quartets = algorithms.quartets_from_splitsystem(system)

API Reference
-------------

**Class**: :class:`phylozoo.core.split.SplitSystem`

**Properties:**

* **elements** - Frozen set of all elements covered by splits.
* **splits** - Frozen set of all splits (read-only).

**Class**: :class:`phylozoo.core.split.WeightedSplitSystem`

**Methods:**

* **get_weight(split)** - Get weight for a specific split. Returns float.

**Algorithms** (``phylozoo.core.split.algorithms``):

* **splitsystem_to_tree(system, check_compatibility=True)** - Convert compatible split 
  system to tree (SemiDirectedPhyNetwork). Returns SemiDirectedPhyNetwork. See 
  :func:`phylozoo.core.split.algorithms.splitsystem_to_tree`.

* **distances_from_splitsystem(system)** - Extract distance matrix from split system. 
  Returns DistanceMatrix. See :func:`phylozoo.core.split.algorithms.distances_from_splitsystem`.

* **quartets_from_splitsystem(system)** - Extract quartet profiles from split system. 
  Returns QuartetProfileSet. See :func:`phylozoo.core.split.algorithms.quartets_from_splitsystem`.

**Classifications** (``phylozoo.core.split.classifications``):

* **is_tree_compatible(system)** - Check if split system is tree-compatible (all splits 
  are pairwise compatible). Returns boolean.

**Conversion Functions:**

* **to_weightedsplitsystem(system, default_weight=1.0)** - Convert SplitSystem to 
  WeightedSplitSystem with default weights. Returns WeightedSplitSystem.

.. tip::
   Use split compatibility to check if a collection of splits can be represented as a 
   tree. If splits are compatible, you can convert them to a tree using 
   ``splitsystem_to_tree``.

.. seealso::
   For single splits, see :doc:`Split <split>`. 
   For extracting splits from networks, see :doc:`Networks (Advanced) <../networks/directed/advanced>`. 
   For I/O operations, see :doc:`I/O <../../io>`.
