Splits
======

Splits represent bipartitions of taxa and are a fundamental way to encode phylogenetic 
relationships :cite:`PhyloZoo2024`. The split module provides classes for working with 
splits and split systems.

Split Basics
------------

A ``Split`` represents a bipartition of a set of elements:

.. code-block:: python

   from phylozoo import Split
   
   # Create a split
   split = Split({"A", "B"}, {"C", "D"})
   
   # Check if trivial (one side has size 1)
   is_trivial = split.is_trivial()
   
   # Check compatibility with another split
   other_split = Split({"A"}, {"B", "C", "D"})
   are_compatible = split.is_compatible(other_split)
   
   # Check if one split is a refinement of another
   is_refinement = split.is_refinement(other_split)
   
   # Check if one split is a subsplit of another
   is_subsplit = split.is_subsplit(other_split)

Split Systems
-------------

A ``SplitSystem`` is a collection of splits where each split covers the complete set 
of elements:

.. code-block:: python

   from phylozoo import SplitSystem
   
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

Complete Workflow Example
--------------------------

Here's a complete example showing how to work with splits:

.. code-block:: python

   from phylozoo import Split, SplitSystem
   from phylozoo.core.split import algorithms
   
   # Create splits representing phylogenetic relationships
   splits = [
       Split({"A", "B"}, {"C", "D", "E"}),  # A,B vs C,D,E
       Split({"A"}, {"B", "C", "D", "E"}),  # A vs others
       Split({"A", "B", "C"}, {"D", "E"}),  # A,B,C vs D,E
   ]
   
   # Create split system
   system = SplitSystem(splits)
   
   # Check if splits are tree-compatible
   from phylozoo.core.split.classifications import is_tree_compatible
   compatible = is_tree_compatible(system)
   
   if compatible:
       # Convert to tree
       tree = algorithms.splitsystem_to_tree(system)
       
       # Extract distance matrix
       distances = algorithms.distances_from_splitsystem(system)
       
       # Extract quartets
       quartets = algorithms.quartets_from_splitsystem(system)

Available Functions
-------------------

**Classes:**

* **Split** - Immutable split class representing a 2-partition of a set. Inherits from 
  Partition. Supports compatibility checking, refinement checking, and subsplit 
  checking. See :class:`phylozoo.core.split.Split` for full API.

* **SplitSystem** - Immutable collection of splits where each split covers the complete 
  set of elements. Validates that all splits cover the same element set. Supports I/O 
  operations. See :class:`phylozoo.core.split.SplitSystem` for full API.

* **WeightedSplitSystem** - Immutable split system with weights associated with each 
  split. Supports weight access and conversion from regular split systems. See 
  :class:`phylozoo.core.split.WeightedSplitSystem` for full API.

**Split Methods:**

* **is_trivial()** - Check if split is trivial (one side has size 1). Returns boolean.
* **is_compatible(other)** - Check compatibility with another split. Two splits are 
  compatible if one side of each split is a subset of one side of the other. Returns 
  boolean.
* **is_refinement(other)** - Check if this split is a refinement of another split. 
  Returns boolean.
* **is_subsplit(other)** - Check if this split is a subsplit of another split. Returns 
  boolean.
* **elements** - Frozen set of all elements in the split.

**SplitSystem Methods:**

* **elements** - Frozen set of all elements covered by splits.
* **splits** - Frozen set of all splits (read-only).

**Utility Functions** (``phylozoo.core.split.base``):

* **is_compatible(split1, split2)** - Check if two splits are compatible. Standalone 
  function version. Returns boolean.

* **is_subsplit(split1, split2)** - Check if split1 is a subsplit of split2. Standalone 
  function version. Returns boolean.

**Algorithms** (``phylozoo.core.split.algorithms``):

* **splitsystem_to_tree(system, check_compatibility=True)** - Convert compatible split 
  system to tree (SemiDirectedPhyNetwork). Uses star tree approach with iterative 
  refinement. Returns SemiDirectedPhyNetwork. See 
  :func:`phylozoo.core.split.algorithms.splitsystem_to_tree`.

* **distances_from_splitsystem(system)** - Extract distance matrix from split system. 
  Computes distances based on split weights. Returns DistanceMatrix. See 
  :func:`phylozoo.core.split.algorithms.distances_from_splitsystem`.

* **quartets_from_splitsystem(system)** - Extract quartet profiles from split system. 
  Groups splits into quartets. Returns QuartetProfileSet. See 
  :func:`phylozoo.core.split.algorithms.quartets_from_splitsystem`.

**Classifications** (``phylozoo.core.split.classifications``):

* **is_tree_compatible(system)** - Check if split system is tree-compatible (all splits 
  are pairwise compatible). Returns boolean.

**Conversion Functions:**

* **to_weightedsplitsystem(system, default_weight=1.0)** - Convert SplitSystem to 
  WeightedSplitSystem with default weights. Returns WeightedSplitSystem.

.. note::
   Split systems can be extracted from networks using derivation functions. See 
   :doc:`Networks (Advanced) <core_networks_advanced>` for details.

.. tip::
   Use split compatibility to check if a collection of splits can be represented as a 
   tree. If splits are compatible, you can convert them to a tree using 
   ``splitsystem_to_tree``.

.. seealso::
   For extracting splits from networks, see :doc:`Networks (Advanced) <core_networks_advanced>`. 
   For I/O operations, see :doc:`I/O <io>`.
