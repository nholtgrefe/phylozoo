Partitions
==========

The :mod:`phylozoo.core.primitives` module provides the :class:`Partition` class, which
represents a mathematical partition of a set into disjoint subsets. Partitions are
fundamental to phylogenetic analysis, serving as the foundation for other classes
in PhyloZoo.

The class can be imported from the core primitives module:

.. code-block:: python

   from phylozoo.core.primitives import *
   # or directly
   from phylozoo.core.primitives import Partition

Working with Partitions
-----------------------

The :class:`phylozoo.core.primitives.Partition` class represents a partition of a set
into disjoint, non-empty subsets that cover the entire set. Partitions are immutable
and provide methods for mathematical operations essential to phylogenetic analysis.

.. note::
   :class: dropdown

   **Implementation details**

   Partitions are designed for mathematical rigor and performance:

   - Internal representation uses frozensets for immutability
   - Elements are automatically sorted for consistent ordering and equality
   - Validation ensures mathematical partition properties (disjoint, covering)
   - Efficient equality and hashing based on canonical representation
   - Memory-efficient storage suitable for large phylogenetic datasets

   For implementation details, see :mod:`src/phylozoo/core/primitives/partition.py`.

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
All subsets must be disjoint and their union must equal the universal set.

Accessing Partition Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Partitions provide several properties for accessing their mathematical components:

.. code-block:: python

   # Access the partition parts
   parts = partition.parts  # frozenset({frozenset({1, 2}), frozenset({3, 4}), frozenset({5})})

   # Get all elements in the partition
   elements = partition.elements  # frozenset({1, 2, 3, 4, 5})

   # Get the number of parts
   num_parts = len(partition)  # 3


Partition Operations
^^^^^^^^^^^^^^^^^^^^

Partitions support fundamental mathematical operations used in phylogenetic analysis:

**Refinement Checking**

The refinement relationship is essential for analyzing hierarchical structures:

.. code-block:: python

   # Create a finer partition (more parts)
   finer = Partition([{1}, {2}, {3, 4}, {5}])

   # Check if finer is a refinement of partition
   is_refinement = partition.is_refinement(finer)  # False
   is_refinement = finer.is_refinement(partition)  # True

A partition P is a refinement of partition Q if every part of P is contained within
some part of Q. Refinement relationships define the hierarchical structure of
phylogenetic trees and networks.

**Mathematical Properties**

Partitions are immutable and hashable, making them suitable for use in sets and
dictionaries. Partitions with the same parts are considered equal regardless of
the order in which parts are specified.

.. code-block:: python

   # Two partitions are equal if they have the same parts
   p1 = Partition([{1, 2}, {3, 4}])
   p2 = Partition([{3, 4}, {1, 2}])  # Same parts, different order
   are_equal = p1 == p2  # True

   # Partitions are hashable (can be used as dict keys or set elements)
   partition_set = {partition, finer}
   partition_dict = {partition: "coarse", finer: "fine"}

**Element Queries**

The :meth:`get_part` method returns the part containing a specific element, and
you can check if a set is one of the partition parts using the `in` operator.

.. code-block:: python

   # Get the part containing a specific element
   part_containing_3 = partition.get_part(3)  # frozenset({3, 4})

   # Check if a set is one of the partition parts
   is_part = {1, 2} in partition  # True
   is_part = {1, 3} in partition  # False

**Size and Counting**

The :meth:`size` method returns the total number of elements across all parts,
while `len()` returns the number of parts in the partition.

.. code-block:: python

   # Number of parts in the partition
   num_parts = len(partition)  # 3

   # Total number of elements across all parts
   total_elements = partition.size()  # 5

**Iteration and Traversal**

Partitions support iteration over all parts, allowing you to process each part
in the partition.

.. code-block:: python

   # Iterate over all parts
   for part in partition:
       print(f"Part: {part}")

   # Convert to list of parts
   parts_list = list(partition)

See Also
--------

- :doc:`API Reference <../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Splits <../splits/overview>` - Split systems based on partitions
- :doc:`Circular Ordering <circular_ordering>` - Circular arrangements for quartet analysis
- :doc:`Networks (Basic) <../networks/basic>` - Network classes using partitions
