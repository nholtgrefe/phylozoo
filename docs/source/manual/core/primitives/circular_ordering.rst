Circular Orderings
==================

Working with Circular Orderings
-------------------------------

Circular orderings represent mathematical arrangements where elements are placed
in a circle, with the first and last elements considered adjacent. PhyloZoo offers two classes:
:class:`CircularOrdering` for individual elements and :class:`CircularSetOrdering`
for grouped elements, both stored in canonical form for efficient mathematical operations.

:class:`CircularSetOrdering` inherits from :class:`Partition`, representing a partition
with an additional circular ordering structure. :class:`CircularOrdering` is a special
case of :class:`CircularSetOrdering` where each element is in its own singleton set.

Creating Circular Orderings
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Circular orderings can be created from any sequence of elements or sets:

.. code-block:: python

   from phylozoo.core.primitives import CircularOrdering, CircularSetOrdering

   # Create circular ordering of individual elements
   co = CircularOrdering([1, 2, 3, 4])

   # Create circular ordering of sets (more general)
   cso = CircularSetOrdering([{1, 2}, {3}, {4, 5}])

   # From different input types
   co_from_tuple = CircularOrdering((1, 2, 3, 4))
   cso_from_list = CircularSetOrdering([[1, 2], [3], [4, 5]])

The constructor automatically converts inputs to the appropriate immutable types
and validates that all elements are present and distinct.

Accessing Circular Ordering Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Circular orderings provide properties to access their structure and elements:

.. code-block:: python

   # Access the element order (CircularOrdering)
   order = co.order  # (1, 2, 3, 4)

   # Access the set order (CircularSetOrdering)
   set_order = cso.setorder  # ({1, 2}, {3}, {4, 5})

   # Get all elements
   elements = co.elements  # frozenset({1, 2, 3, 4})
   set_elements = cso.elements  # frozenset({1, 2, 3, 4, 5})

   # Get number of components
   num_elements = len(co)  # 4
   num_sets = len(cso)     # 3

Canonical Representations
^^^^^^^^^^^^^^^^^^^^^^^^^

Circular orderings use canonical representations for efficient equality comparisons.
Two orderings are considered equal if they represent the same circular arrangement,
regardless of starting position or direction:

.. code-block:: python

   # Different representations of the same circular ordering
   co1 = CircularOrdering([1, 2, 3, 4])
   co2 = CircularOrdering([2, 3, 4, 1])  # Rotated
   co3 = CircularOrdering([4, 3, 2, 1])  # Reversed

   # All represent the same circular arrangement
   co1 == co2  # True
   co1 == co3  # True
   co2 == co3  # True

   # Hashable for use in sets and dicts
   ordering_set = {co1, co2, co3}  # Contains only one unique ordering

This canonical form enables efficient mathematical operations and eliminates
duplicates in collections of circular orderings.

Circular Ordering Operations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary operation on circular orderings is neighbor checking, which determines
whether two elements or sets are adjacent in the circular arrangement:

.. code-block:: python

   # Check element neighbors (CircularOrdering)
   adjacent = co.are_neighbors(1, 2)    # True (consecutive)
   wrapped = co.are_neighbors(1, 4)     # True (wraps around)
   not_adjacent = co.are_neighbors(1, 3) # False (not neighbors)

   # Check set neighbors (CircularSetOrdering)
   set_adjacent = cso.are_neighbors({1, 2}, {3})        # True
   set_wrapped = cso.are_neighbors({1, 2}, {4, 5})      # True (wraps around)
   set_not_adjacent = cso.are_neighbors({1, 2}, {4, 5}) # False in this case

Neighbor relationships are symmetric and consider the circular topology where
the first and last elements are adjacent.


See Also
--------

- :doc:`API Reference <../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Partition <partition>` - Set partitions for hierarchical analysis