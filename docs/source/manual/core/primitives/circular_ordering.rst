Circular Ordering
=================

Circular orderings represent elements arranged in a circle. PhyloZoo provides two 
classes: ``CircularOrdering`` for elements and ``CircularSetOrdering`` for sets.

CircularOrdering
----------------

A ``CircularOrdering`` represents a circular arrangement of individual elements:

.. code-block:: python

   from phylozoo.core.primitives import CircularOrdering
   
   # Create circular ordering of elements
   co = CircularOrdering([1, 2, 3, 4])
   
   # Check if elements are neighbors
   are_neighbors = co.are_neighbors(1, 2)  # True
   are_neighbors = co.are_neighbors(1, 4)  # True (circular, wraps around)
   
   # Access order
   order = co.order  # Tuple of elements in order

CircularSetOrdering
-------------------

A ``CircularSetOrdering`` represents a circular arrangement of sets (more general 
than CircularOrdering):

.. code-block:: python

   from phylozoo.core.primitives import CircularSetOrdering
   
   # Create circular ordering of sets
   cso = CircularSetOrdering([{1, 2}, {3}, {4, 5}])
   
   # Check if sets are neighbors
   are_neighbors = cso.are_neighbors({1, 2}, {3})  # True
   
   # Access set order
   setorder = cso.setorder  # Tuple of sets in order

Canonical Form
--------------

Circular orderings are stored in canonical form, making equality comparisons 
efficient. Two orderings are equal if they are the same up to cyclic permutation 
or reversal:

.. code-block:: python

   co1 = CircularOrdering([1, 2, 3, 4])
   co2 = CircularOrdering([2, 3, 4, 1])  # Cyclic permutation
   co3 = CircularOrdering([4, 3, 2, 1])  # Reversal
   
   co1 == co2  # True
   co1 == co3  # True

API Reference
-------------

**Class**: :class:`phylozoo.core.primitives.CircularOrdering`

**Properties:**

* **order** - Tuple of elements in circular order

**Methods:**

* **are_neighbors(a, b)** - Check if two elements are neighbors in the circular ordering. Returns boolean.

**Class**: :class:`phylozoo.core.primitives.CircularSetOrdering`

**Properties:**

* **setorder** - Tuple of sets in circular order

**Methods:**

* **are_neighbors(set1, set2)** - Check if two sets are neighbors in the circular ordering. Returns boolean.

.. tip::
   Circular orderings are particularly useful for quartet-based methods and distance 
   matrix classifications (e.g., Kalmanson matrices).

.. seealso::
   For quartets (which use circular orderings), see :doc:`Quartets <../quartets/overview>`.
