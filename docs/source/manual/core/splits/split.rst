Splits
======

Working with Splits
--------------------

The :mod:`phylozoo.core.split` module provides the :class:`~phylozoo.core.split.base.Split` class, which represents
a bipartition of a set of elements into two complementary subsets. Splits are fundamental
to phylogenetic analysis, providing the mathematical foundation for representing
evolutionary relationships in trees and networks. The split class inherits from 
the :class:`~phylozoo.core.primitives.partition.Partition` class, which provides additional 
mathematical operations.


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
is properly formed. It raises a :class:`~phylozoo.utils.exceptions.general.PhyloZooValueError` if the input is invalid.

Accessing Split Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

Splits provide access to their components:

.. code-block:: python

   # Access all elements in the split
   elements = split.elements  # frozenset({"A", "B", "C", "D"})

   # Access individual parts (inherited from Partition)
   parts = split.parts        # tuple of frozensets
   part1, part2 = split.parts # (frozenset({"A", "B"}), frozenset({"C", "D"}))

Split Operations
--------------------

Splits support fundamental phylogenetic operations.

**Equality Checking**

The `==` operator (which is inherited from :class:`~phylozoo.core.primitives.partition.Partition`) checks whether two splits are equal. 
Two splits are equal if they have the same parts, regardless of the order of the parts.

.. code-block:: python

   # Check if two splits are equal
   split1 = Split({"A", "B"}, {"C", "D"})
   split2 = Split({"C", "D"}, {"A", "B"})
   are_equal = split1 == split2  # True

**Triviality Checking**

The :meth:`~phylozoo.core.split.base.Split.is_trivial` method checks whether a split is trivial, meaning one side
contains only a single element.

.. code-block:: python

   # Check if split is trivial (one side has only one element)
   trivial_split = Split({"A"}, {"B", "C", "D"})
   is_trivial = trivial_split.is_trivial()  # True

   normal_split = Split({"A", "B"}, {"C", "D"})
   is_trivial = normal_split.is_trivial()  # False

**Compatibility Analysis**

The :func:`~phylozoo.core.split.classifications.is_compatible` function checks whether two splits are
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

**Subsplit Relationships**

The :func:`~phylozoo.core.split.classifications.is_subsplit` function checks whether one split is a
subsplit of another. A split :math:`S_1` is a subsplit of :math:`S_2` if one of its 
sides is a subset of one side of the other split, and the other side of this split 
is a subset of the other side of the other split. For example, :math:`\{1, 2\} | \{5, 6\}` is a 
subsplit of :math:`\{1, 2, 3\} | \{4, 5, 6\}` because :math:`\{1, 2\} \subseteq \{1, 2, 3\}` and :math:`\{5, 6\} \subseteq \{4, 5, 6\}`.

.. code-block:: python

   from phylozoo.core.split import is_subsplit

   # Check subsplit relationships
   split = Split({"1", "2", "6"}, {"3", "4", "5"})
   subsplit = Split({"1", "2"}, {"4", "5"})

   is_subsplit = is_subsplit(subsplit, split)  # True

The :func:`~phylozoo.core.split.algorithms.induced_quartetsplits` function returns a set of all subsplits of size 4 of the split.

.. code-block:: python

   from phylozoo.core.split import induced_quartetsplits

   # Return all quartets splits
   q_splits = induced_quartetsplits(split)

See Also
--------

- :doc:`API Reference <../../../../api/core/splits>` - Complete function signatures and detailed examples
- :doc:`Split Systems <split_system>` - Collections of splits for phylogenetic analysis
- :doc:`Partitions <../primitives/partition>` - Base partition class
- :doc:`Quartets <../quartets/overview>` - Quartet-based phylogenetic representations
