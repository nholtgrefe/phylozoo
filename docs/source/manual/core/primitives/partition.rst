Partition
==========

A ``Partition`` represents a partition of a set into disjoint subsets. It is used as 
a base class for ``Split`` and provides methods for checking refinement and other 
partition operations.

Creating Partitions
-------------------

.. code-block:: python

   from phylozoo.core.primitives import Partition
   
   # Create a partition
   partition = Partition([{1, 2}, {3, 4}, {5}])
   
   # Access parts
   parts = partition.parts  # Frozen set of frozensets
   elements = partition.elements  # All elements in the partition

Partition Operations
--------------------

.. code-block:: python

   # Check if valid partition
   is_valid = partition.is_partition()  # True
   
   # Check refinement
   other = Partition([{1}, {2}, {3, 4}, {5}])
   is_refinement = partition.is_refinement(other)  # True
   
   # Get size (number of parts)
   size = len(partition)  # 3

API Reference
-------------

**Class**: :class:`phylozoo.core.primitives.Partition`

**Properties:**

* **parts** - Frozen set of frozensets representing the partition parts
* **elements** - Frozen set of all elements in the partition

**Methods:**

* **is_partition()** - Check if this is a valid partition (disjoint subsets covering all elements). Returns boolean.
* **is_refinement(other)** - Check if this partition is a refinement of another partition. Returns boolean.

**Size:**

* **len(partition)** - Number of parts in the partition

.. seealso::
   For splits (which inherit from Partition), see :doc:`Split <../splits/split>`.
