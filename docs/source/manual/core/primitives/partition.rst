Partitions
==========

Working with Partitions
-----------------------

The :class:`~phylozoo.core.primitives.partition.Partition` class represents a partition of a set
into disjoint, non-empty subsets that cover the entire set. Partitions are immutable
and provide methods for mathematical operations essential to phylogenetic analysis.

Creating Partitions
^^^^^^^^^^^^^^^^^^^

Partitions can be created from any iterable of iterables representing the disjoint subsets:

.. code-block:: python

   from phylozoo.core.primitives import Partition

   # Create a partition from sets
   partition = Partition([{1, 2}, {3, 4}, {5}])

   # Create from lists (automatically converted to sets)
   partition2 = Partition([[1, 2], [3, 4], [5]])

   # Empty partition (edge case)
   empty_partition = Partition([])

The constructor validates that the input represents a valid mathematical partition.
All subsets must be disjoint.

Accessing Partition Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Partitions provide access to their parts and the set of all elements:

.. code-block:: python

   # Access the partition parts
   parts = partition.parts  # frozenset({frozenset({1, 2}), frozenset({3, 4}), frozenset({5})})

   # Get all elements in the partition
   elements = partition.elements  # frozenset({1, 2, 3, 4, 5})

**Equality and hashing**

Partitions are immutable and hashable. Two partitions are considered equal if they
have the same parts, regardless of the order in which parts were specified:

.. code-block:: python

   # Two partitions are equal if they have the same parts
   p1 = Partition([{1, 2}, {3, 4}])
   p2 = Partition([{3, 4}, {1, 2}])  # Same parts, different order
   are_equal = p1 == p2  # True

   # Partitions are hashable (can be used as dict keys or set elements)
   partition_set = {partition, finer}
   partition_dict = {partition: "coarse", finer: "fine"}

**Element queries**

The :meth:`~phylozoo.core.primitives.partition.Partition.get_part` method returns the part
containing a specific element. You can check if a set is one of the partition parts
using the ``in`` operator:

.. code-block:: python

   # Get the part containing a specific element
   part_containing_3 = partition.get_part(3)  # frozenset({3, 4})

   # Check if a set is one of the partition parts
   is_part = {1, 2} in partition  # True
   is_part = {1, 3} in partition  # False

**Size and length**

:meth:`~phylozoo.core.primitives.partition.Partition.size` returns the total number of
elements; :func:`len` returns the number of parts:

.. code-block:: python

   # Number of parts in the partition
   num_parts = len(partition)  # 3

   # Total number of elements across all parts
   total_elements = partition.size()  # 5

**Iteration**

You can iterate over the parts of a partition (in canonical order):

.. code-block:: python

   for part in partition:
       print(part)  # Each part is a frozenset

Partition Operations
--------------------

Partitions support fundamental mathematical operations used in phylogenetic analysis.

**Refinement checking**

A partition :math:`P` is a refinement of partition :math:`Q` if every part of :math:`P` is contained within
some part of :math:`Q`. Refinement relations can be checked using the :meth:`~phylozoo.core.primitives.partition.Partition.is_refinement` method:

.. code-block:: python

   # Create a finer partition (more parts)
   finer = Partition([{1}, {2}, {3, 4}, {5}])

   # Check if finer is a refinement of partition
   is_refinement = partition.is_refinement(finer)  # False
   is_refinement = finer.is_refinement(partition)  # True

**Subpartitions**

The :meth:`~phylozoo.core.primitives.partition.Partition.subpartitions` method yields all subpartitions of a given
size. A subpartition is formed by selecting that many parts from the current partition.

.. code-block:: python

   # All subpartitions with 2 parts from a 5-part partition
   partition = Partition([{1}, {2}, {3}, {4}, {5}])
   subparts = list(partition.subpartitions(size=2))
   len(subparts)  # C(5, 2) = 10

**Representative partitions**

The :meth:`~phylozoo.core.primitives.partition.Partition.representative_partitions` method yields all partitions
that have exactly one element chosen from each part (singleton representatives). Useful for
enumerating one representative per block.

.. code-block:: python

   partition = Partition([{1, 2}, {3, 4}])
   reps = list(partition.representative_partitions())
   len(reps)  # 2 * 2 = 4 (one element from first part, one from second)
   reps[0]    # e.g. Partition([{1}, {3}])

See Also
--------

- :doc:`API Reference <../../../../api/core/primitives/index>` - Complete function signatures and detailed examples
- :doc:`Splits <../splits/overview>` - A specific partition with exactly two parts
- :doc:`Circular Ordering <circular_ordering>` - A partition with a circular ordering structure
