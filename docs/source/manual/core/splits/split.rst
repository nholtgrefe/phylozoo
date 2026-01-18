Split
=====

A ``Split`` represents a bipartition of a set of elements. It inherits from ``Partition`` 
and provides methods for checking compatibility, refinement, and subsplit relationships.

Creating Splits
---------------

.. code-block:: python

   from phylozoo import Split
   
   # Create a split
   split = Split({"A", "B"}, {"C", "D"})

Split Operations
----------------

.. code-block:: python

   # Check if trivial (one side has size 1)
   is_trivial = split.is_trivial()
   
   # Check compatibility with another split
   other_split = Split({"A"}, {"B", "C", "D"})
   are_compatible = split.is_compatible(other_split)
   
   # Check if one split is a refinement of another
   is_refinement = split.is_refinement(other_split)
   
   # Check if one split is a subsplit of another
   is_subsplit = split.is_subsplit(other_split)
   
   # Access elements
   elements = split.elements  # Frozen set of all elements

API Reference
-------------

**Class**: :class:`phylozoo.core.split.Split`

**Methods:**

* **is_trivial()** - Check if split is trivial (one side has size 1). Returns boolean.
* **is_compatible(other)** - Check compatibility with another split. Two splits are 
  compatible if one side of each split is a subset of one side of the other. Returns boolean.
* **is_refinement(other)** - Check if this split is a refinement of another split. 
  Returns boolean.
* **is_subsplit(other)** - Check if this split is a subsplit of another split. Returns boolean.

**Properties:**

* **elements** - Frozen set of all elements in the split.

**Utility Functions** (``phylozoo.core.split.base``):

* **is_compatible(split1, split2)** - Check if two splits are compatible. Standalone 
  function version. Returns boolean.

* **is_subsplit(split1, split2)** - Check if split1 is a subsplit of split2. Standalone 
  function version. Returns boolean.

.. seealso::
   For collections of splits, see :doc:`SplitSystem <split_system>`.
