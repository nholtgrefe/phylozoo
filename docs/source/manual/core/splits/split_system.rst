Split Systems
=============

Working with Split Systems
--------------------------

The :mod:`phylozoo.core.split` module provides the :class:`~phylozoo.core.split.splitsystem.SplitSystem` and :class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem` classes,
which represent collections of splits covering a complete set of elements.
Split systems provide a comprehensive way to represent phylogenetic relationships
across multiple taxa.

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

The constructor validates that all splits cover the same set of elements and that there are no duplicate splits. It raises a :class:`~phylozoo.utils.exceptions.general.PhyloZooValueError` if the input is invalid.

Creating Weighted Split Systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Weighted split systems extend regular split systems with weights, enabling, e.g., 
phylogenetic analyses where splits have different levels of support.

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

Accessing Properties
^^^^^^^^^^^^^^^^^^^^^

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

Weighted split systems provide additional methods to access split weights:

.. code-block:: python

   # Get weight for a specific split
   weight = weighted_system.get_weight(splits[0])

   # Get all weights as a dictionary
   all_weights = weighted_system.weights

   # Get total weight sum
   total = weighted_system.total_weight

File Input/Output
^^^^^^^^^^^^^^^^^

Split systems support reading and writing in NEXUS format:

- **NEXUS** (default): Standard phylogenetic data format for split systems — see :doc:`NEXUS format <../../utils/io/formats/nexus>`

.. code-block:: python

   # Load from file (auto-detects format by extension)
   split_system = SplitSystem.load("splits.nexus")

   # Load with explicit format
   split_system = SplitSystem.load("splits.nex", format="nexus")

   # Save to file
   split_system.save("output.nexus")

Weighted split systems also support NEXUS format with weights:

.. code-block:: python

   # Load weighted split system
   weighted_system = WeightedSplitSystem.load("weighted_splits.nexus")

   # Save weighted split system
   weighted_system.save("output.nexus")

.. seealso::
   The `SplitSystem` and `WeightedSplitSystem` classes use the :class:`~phylozoo.utils.io.mixin.IOMixin` interface, providing
   consistent file handling across PhyloZoo classes. For details on the I/O system,
   see the :doc:`I/O manual <../../utils/io/overview>`.

Classifications
---------------

The split module provides functions to classify split systems based on mathematical properties. 
These classifications are essential for validating that split systems meet the requirements of specific algorithms.

**Trivial Splits**

The :func:`~phylozoo.core.split.classifications.has_all_trivial_splits` function checks
whether the split system contains all trivial splits. See :meth:`~phylozoo.core.split.base.Split.is_trivial` for the definition of trivial splits. 
For a split system with :math:`n` elements, there should be exactly :math:`n` trivial splits.

.. code-block:: python

   from phylozoo.core.split.classifications import has_all_trivial_splits

   # Check if all trivial splits are present
   has_all = has_all_trivial_splits(split_system)

**Pairwise Compatibility**

The :func:`~phylozoo.core.split.classifications.is_pairwise_compatible` function checks
whether all pairs of splits in the system are compatible with each other. See 
:func:`~phylozoo.core.split.classifications.is_compatible` for the definition of compatibility.

.. code-block:: python

   from phylozoo.core.split.classifications import is_pairwise_compatible

   # Check if all pairs of splits are compatible
   compatible = is_pairwise_compatible(split_system)

A split system is pairwise compatible if every pair of splits in the system is
compatible with each other.

**Tree Compatibility**

The :func:`~phylozoo.core.split.classifications.is_tree_compatible` function checks
whether a split system is compatible with a tree. A split system is tree-compatible if:
1. All pairs of splits are compatible (pairwise compatible)
2. All trivial splits are present in the system

.. code-block:: python

   from phylozoo.core.split.classifications import is_tree_compatible

   # Check if split system can form a phylogenetic tree
   compatible = is_tree_compatible(split_system)

   if compatible:
       print("Split system is tree-compatible")
   else:
       print("Split system contains incompatible splits")

Algorithms
----------

Split systems support fundamental phylogenetic operations for analysis and conversion.

Tree Reconstruction
^^^^^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.split.algorithms.tree_from_splitsystem` function converts a
compatible split system into an unrooted phylogenetic tree represented as a `~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`.
The algorithm builds a tree that displays all splits in the system.

.. code-block:: python

   from phylozoo.core.split.algorithms import tree_from_splitsystem

   # Convert compatible split system to phylogenetic tree
   if is_tree_compatible(split_system):
       tree = tree_from_splitsystem(split_system)

The function requires the split system to be tree-compatible. If compatibility is not
guaranteed, the function will raise an error unless `check_compatibility=False` is
specified.

Distance Matrix Extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.split.algorithms.distances_from_splitsystem` function computes
a distance matrix from a split system. The distance between two taxa is
the sum of weights of all splits that separate them.

.. code-block:: python

   from phylozoo.core.split.algorithms import distances_from_splitsystem

   # Extract evolutionary distance matrix from split system
   distance_matrix = distances_from_splitsystem(split_system)

For weighted split systems, the split weights are used in the distance computation.
For regular split systems, each split contributes weight 1.0 to the distances.

Quartet Profile Extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.split.algorithms.quartets_from_splitsystem` function extracts
quartet profiles (see :doc:`Quartets <../quartets/overview>`) from a split system. For each split, the function generates all
quartets induced by that split, then groups them by four-taxon sets into profiles.

.. code-block:: python

   from phylozoo.core.split.algorithms import quartets_from_splitsystem

   # Extract quartet profiles from split system
   quartet_profiles = quartets_from_splitsystem(split_system)

See Also
--------

- :doc:`API Reference <../../../api/core/splits>` - Complete function signatures and detailed examples
- :doc:`Splits <split>` - Individual bipartitions
- :doc:`Quartets <../quartets/overview>` - Quartet-based representations
- :doc:`Distance Matrices <../distance>` - Distance matrix computations
