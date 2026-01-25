Splits
======

The :mod:`phylozoo.core.split` module provides the :class:`Split` class, which represents
a bipartition of a set of taxa into two complementary subsets. Splits are fundamental
to phylogenetic analysis, providing the mathematical foundation for representing
evolutionary relationships in trees and networks.

All classes and functions on this page can be imported from the core split module:

.. code-block:: python

   from phylozoo.core.split import *
   # or directly
   from phylozoo.core.split import Split

Working with Splits
--------------------

Splits represent the fundamental phylogenetic concept of dividing taxa into two
complementary groups. They form the mathematical basis for phylogenetic trees and
networks, where each internal edge corresponds to a split of the taxa.

.. note::
   :class: dropdown

   **Implementation details**

   Split representation is optimized for phylogenetic operations:

   - Canonical ordering ensures consistent equality and hashing
   - Efficient compatibility checking algorithms
   - Immutable design ensuring mathematical consistency
   - Inheritance from Partition provides additional set operations

   For implementation details, see :mod:`src/phylozoo/core/split/split.py`.

Creating Splits
^^^^^^^^^^^^^^^^

Splits can be created from any two complementary subsets:

.. code-block:: python

   from phylozoo.core.split import Split

   # Create a split dividing taxa into two groups
   split = Split({"A", "B"}, {"C", "D"})

   # Create from different iterable types
   split2 = Split(["A", "B"], ["C", "D"])

The constructor automatically converts inputs to frozensets and ensures the split
is properly formed.

Accessing Split Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

Splits provide access to their mathematical components:

.. code-block:: python

   # Access all elements in the split
   elements = split.elements  # frozenset({"A", "B", "C", "D"})

   # Access individual parts (inherited from Partition)
   parts = split.parts        # tuple of frozensets
   part1, part2 = split.parts # (frozenset({"A", "B"}), frozenset({"C", "D"}))

Split Operations
^^^^^^^^^^^^^^^^

Splits support fundamental phylogenetic operations.

**Triviality Checking**

The :meth:`is_trivial` method checks whether a split is trivial, meaning one side
contains only a single element. Trivial splits correspond to external edges in
phylogenetic trees.

.. code-block:: python

   # Check if split is trivial (one side has only one element)
   trivial_split = Split({"A"}, {"B", "C", "D"})
   is_trivial = trivial_split.is_trivial()  # True

   normal_split = Split({"A", "B"}, {"C", "D"})
   is_trivial = normal_split.is_trivial()  # False

**Compatibility Analysis**

The :func:`phylozoo.core.split.is_compatible` function checks whether two splits are
compatible, meaning they can coexist in the same phylogenetic tree. Two splits are
compatible if one side of each split is a subset of one side of the other split.

.. code-block:: python

   from phylozoo.core.split import is_compatible

   # Check compatibility between splits
   split1 = Split({"A", "B"}, {"C", "D"})
   split2 = Split({"A"}, {"B", "C", "D"})
   compatible = is_compatible(split1, split2)  # True

   split3 = Split({"A", "C"}, {"B", "D"})
   compatible = is_compatible(split1, split3)  # False (incompatible)

**Refinement Relationships**

The :meth:`is_refinement` method checks whether one split is a refinement of another,
meaning it represents a finer division of the taxa. A split S1 is a refinement of
S2 if one side of S1 is a proper subset of one side of S2.

.. code-block:: python

   # Check refinement relationships
   coarse = Split({"A", "B"}, {"C", "D", "E"})
   fine = Split({"A"}, {"B", "C", "D", "E"})

   is_refinement = fine.is_refinement(coarse)  # True
   is_refinement = coarse.is_refinement(fine)  # False

**Subsplit Relationships**

The :func:`phylozoo.core.split.is_subsplit` function checks whether one split is a
subsplit of another. A split S1 is a subsplit of S2 if one side of S1 is a subset
of one side of S2. This relationship is used in hierarchical phylogenetic analysis.

.. code-block:: python

   from phylozoo.core.split import is_subsplit

   # Check subsplit relationships
   split = Split({"A", "B"}, {"C", "D"})
   subsplit = Split({"A"}, {"B", "C", "D"})

   is_subsplit = is_subsplit(subsplit, split)  # True

**Mathematical Properties**

Splits are immutable and hashable, making them suitable for use in sets and
dictionaries. Splits with the same topology are considered equal.

.. code-block:: python

   # Splits are immutable and hashable
   split_copy = Split({"A", "B"}, {"C", "D"})
   are_equal = split == split_copy  # True

   # Can be used as dictionary keys or set elements
   split_dict = {split: "data"}
   split_set = {split, split_copy}  # Contains one element

See Also
--------

- :doc:`API Reference <../../../api/core/split>` - Complete function signatures and detailed examples
- :doc:`Split Systems <split_system>` - Collections of splits for phylogenetic analysis
- :doc:`Partitions <../primitives/partition>` - Base partition class
- :doc:`Quartets <../quartets/overview>` - Quartet-based phylogenetic representations
