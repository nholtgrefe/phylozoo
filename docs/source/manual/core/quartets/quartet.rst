Quartet
=======

A ``Quartet`` represents an unrooted tree on four taxa. There are three possible 
quartet topologies, or a star quartet (unresolved).

Creating Quartets
-----------------

.. code-block:: python

   from phylozoo import Quartet
   from phylozoo.core.split import Split
   
   # Create resolved quartets
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))  # AB|CD
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))  # AC|BD
   q3 = Quartet(Split({"A", "D"}, {"B", "C"}))  # AD|BC
   
   # Create star quartet (unresolved)
   star = Quartet({"A", "B", "C", "D"})

Quartet Properties
------------------

.. code-block:: python

   # Check if resolved (has a split)
   is_resolved = q1.is_resolved()
   
   # Check if star (unresolved)
   is_star = star.is_star()
   
   # Check compatibility with another quartet
   are_compatible = q1.is_compatible(q2)
   
   # Access taxa
   taxa = q1.taxa  # Frozen set of 4 taxon labels

API Reference
-------------

**Class**: :class:`phylozoo.core.quartet.Quartet`

**Methods:**

* **is_resolved()** - Check if quartet is resolved (has a split). Returns boolean.
* **is_star()** - Check if quartet is unresolved (star tree). Returns boolean.
* **is_compatible(other)** - Check compatibility with another quartet. Returns boolean.

**Properties:**

* **taxa** - Frozen set of 4 taxon labels.

.. seealso::
   For quartet profiles, see :doc:`QuartetProfile <quartet_profile>`. 
   For collections of quartets, see :doc:`QuartetProfileSet <quartet_profile_set>`.
