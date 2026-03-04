Circular Orderings
==================

Working with Circular Orderings
-------------------------------

Circular orderings represent mathematical arrangements where elements are placed
in a circle, with the first and last elements considered adjacent. PhyloZoo offers two classes:
:class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering` for individual elements
and :class:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering` for grouped
elements, both stored in canonical form for efficient mathematical operations.

:class:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering` inherits from
:class:`~phylozoo.core.primitives.partition.Partition`, representing a partition with an
additional circular ordering structure. :class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering`
is a special case of :class:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering`
where each element is in its own singleton set.

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

**Basic properties**

Circular orderings provide access to their order, elements, and length:

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

**Canonical representations**

Circular orderings use canonical representations for equality. Two orderings are
equal if they represent the same circular arrangement, regardless of starting
position or direction:

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

Circular Ordering Operations
-----------------------------

**Neighbor checking**

The :meth:`~phylozoo.core.primitives.circular_ordering.CircularOrdering.are_neighbors`
method determines whether two elements or sets are adjacent in the circular
arrangement (neighbor relationships are symmetric and wrap around):

.. code-block:: python

   # Check element neighbors (CircularOrdering)
   adjacent = co.are_neighbors(1, 2)    # True (consecutive)
   wrapped = co.are_neighbors(1, 4)     # True (wraps around)
   not_adjacent = co.are_neighbors(1, 3)  # False (not neighbors)

   # Check set neighbors (CircularSetOrdering)
   set_adjacent = cso.are_neighbors({1, 2}, {3})        # True
   set_wrapped = cso.are_neighbors({1, 2}, {4, 5})      # True (wraps around)
   set_not_adjacent = cso.are_neighbors({1, 2}, {4, 5})  # False in this case

**Suborderings**

The :meth:`~phylozoo.core.primitives.circular_ordering.CircularOrdering.suborderings` (respectively, :meth:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering.suborderings`) method yields all
circular (set) orderings formed by selecting a subset of elements (respectively, sets) of a given size:

.. code-block:: python

   # All circular orderings of 2 elements from 5
   co = CircularOrdering([1, 2, 3, 4, 5])
   subs = list(co.suborderings(size=2))
   len(subs)  # C(5, 2) = 10

   # All circular set orderings of 2 sets from 5
   cso = CircularSetOrdering([{1}, {2}, {3}, {4}, {5}])
   set_subs = list(cso.suborderings(size=2))

**Representative orderings**

The :meth:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering.representative_orderings` method yields all
:class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering` instances with exactly one element chosen from each set in the circular set ordering.

.. code-block:: python

   cso = CircularSetOrdering([{1, 2}, {3}])
   reps = list(cso.representative_orderings())
   len(reps)  # 2 (two choices from first set, one from second)

**Converting between CircularOrdering and CircularSetOrdering**

The :meth:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering.are_singletons` method checks if all sets in a :class:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering` are singletons.  If this is the case, the :meth:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering.to_circular_ordering` method converts the :class:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering` to a :class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering`.

.. code-block:: python

   cso = CircularSetOrdering([{1}, {2}, {3}])
   singletons = cso.are_singletons() # True
   co = cso.to_circular_ordering()
 
A :class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering` can be converted to a :class:`~phylozoo.core.primitives.circular_ordering.CircularSetOrdering` using the :meth:`~phylozoo.core.primitives.circular_ordering.CircularOrdering.to_circular_setordering` method.

.. code-block:: python

   co = CircularOrdering([1, 2, 3])
   cso = co.to_circular_setordering()



See Also
--------

- :doc:`API Reference <../../../../api/core/primitives/index>` - Complete function signatures and detailed examples
- :doc:`Partition <partition>` - Set partitions forming the basis of the circular ordering