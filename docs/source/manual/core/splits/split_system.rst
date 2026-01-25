Split Systems
============

The :mod:`phylozoo.core.split` module provides the :class:`SplitSystem` class, which
represents a collection of splits covering a complete set of taxa. Split systems are
fundamental to phylogenetic analysis, encoding complete evolutionary relationships that
can represent trees, networks, or partial phylogenetic information.

All classes and functions on this page can be imported from the core split module:

.. code-block:: python

   from phylozoo.core.split import *
   # or directly
   from phylozoo.core.split import SplitSystem, WeightedSplitSystem

Working with Split Systems
--------------------------

Split systems provide a comprehensive way to represent phylogenetic relationships
across multiple taxa. They validate that all splits cover the same element set and
support various operations for phylogenetic analysis and conversion.

.. note::
   :class: dropdown

   **Implementation details**

   Split systems are optimized for large-scale phylogenetic analysis:

   - Efficient validation of split compatibility and coverage
   - Immutable design ensuring data integrity
   - Support for weighted splits in statistical analyses
   - Integration with tree reconstruction and network inference algorithms

   For implementation details, see :mod:`src/phylozoo/core/split/split_system.py`.

Creating Split Systems
^^^^^^^^^^^^^^^^^^^^^^

Split systems can be created from collections of splits:

.. code-block:: python

   from phylozoo.core.split import SplitSystem, Split

   # Create individual splits
   splits = [
       Split({"A", "B"}, {"C", "D"}),
       Split({"A"}, {"B", "C", "D"}),
       Split({"A", "B", "C"}, {"D"})
   ]

   # Create split system
   split_system = SplitSystem(splits)

The constructor validates that all splits cover the same set of elements.

Accessing Split System Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Split systems provide access to their structure and components:

.. code-block:: python

   # Basic properties
   elements = split_system.elements  # frozenset of all taxa
   num_splits = len(split_system)    # number of splits

   # Access splits
   all_splits = split_system.splits  # frozenset of Split objects

   # Iterate over splits
   for split in split_system:
       print(f"Split: {split}")

   # Check membership
   has_split = some_split in split_system

Split System Operations
-----------------------

Split systems support fundamental phylogenetic operations for analysis and conversion.

Tree Compatibility Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.split.classifications.is_tree_compatible` function checks
whether a split system can be represented as a phylogenetic tree. A split system is
tree-compatible if all pairs of splits are compatible and all trivial splits are present.

.. code-block:: python

   from phylozoo.core.split.classifications import is_tree_compatible

   # Check if split system can form a phylogenetic tree
   compatible = is_tree_compatible(split_system)

   if compatible:
       print("Split system is tree-compatible")
   else:
       print("Split system contains incompatible splits")

Tree Reconstruction
^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.split.algorithms.splitsystem_to_tree` function converts a
compatible split system into a phylogenetic tree represented as a SemiDirectedPhyNetwork.
The algorithm builds a tree that displays all splits in the system.

.. code-block:: python

   from phylozoo.core.split.algorithms import splitsystem_to_tree

   # Convert compatible split system to phylogenetic tree
   if is_tree_compatible(split_system):
       tree = splitsystem_to_tree(split_system)

The function requires the split system to be tree-compatible. If compatibility is not
guaranteed, the function will raise an error unless `check_compatibility=False` is
specified.

Distance Matrix Extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.split.algorithms.distances_from_splitsystem` function computes
an evolutionary distance matrix from a split system. The distance between two taxa is
the sum of weights of all splits that separate them.

.. code-block:: python

   from phylozoo.core.split.algorithms import distances_from_splitsystem

   # Extract evolutionary distance matrix from split system
   distance_matrix = distances_from_splitsystem(split_system)

For weighted split systems, the split weights are used in the distance computation.
For regular split systems, each split contributes weight 1.0 to the distances.

Quartet Profile Extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.split.algorithms.quartets_from_splitsystem` function extracts
quartet profiles from a split system. For each split, the function generates all
quartets induced by that split, then groups them by four-taxon sets into profiles.

.. code-block:: python

   from phylozoo.core.split.algorithms import quartets_from_splitsystem

   # Extract quartet profiles from split system
   quartet_profiles = quartets_from_splitsystem(split_system)

This operation is useful for converting split-based representations to quartet-based
analyses, enabling compatibility with quartet inference algorithms.

Weighted Split Systems
----------------------

Weighted split systems extend regular split systems with confidence scores, enabling
statistical phylogenetic analyses where splits have different levels of support.

Creating Weighted Split Systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Weighted split systems can be created from split-weight pairs or converted from
regular split systems:

.. code-block:: python

   from phylozoo.core.split import WeightedSplitSystem, to_weightedsplitsystem

   # Create weighted split system from list of (split, weight) tuples
   weighted_system = WeightedSplitSystem([
       (splits[0], 0.9),
       (splits[1], 0.8),
       (splits[2], 0.7)
   ])

   # Convert regular split system to weighted
   weighted = to_weightedsplitsystem(split_system, default_weight=1.0)

Accessing Weights
^^^^^^^^^^^^^^^^^

Weighted split systems provide methods to access split weights:

.. code-block:: python

   # Get weight for a specific split
   weight = weighted_system.get_weight(splits[0])

   # Get all weights as a dictionary
   all_weights = weighted_system.weights

   # Get total weight sum
   total = weighted_system.total_weight

See Also
--------

- :doc:`API Reference <../../../api/core/split>` - Complete function signatures and detailed examples
- :doc:`Splits <split>` - Individual bipartitions
- :doc:`Quartets <../quartets/overview>` - Quartet-based representations
- :doc:`Distance Matrices <../distance>` - Distance matrix computations