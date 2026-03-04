Quartets
========

Working with Quartets
----------------------

The :mod:`phylozoo.core.quartet` module provides the :class:`~phylozoo.core.quartet.base.Quartet` class, which
represents an unrooted tree topology on four taxa. Quartets are the fundamental units
of phylogenetic relationships in quartet-based methods, with three possible resolved
topologies or an unresolved star quartet.

Creating Quartets
^^^^^^^^^^^^^^^^^

Quartets can be created from splits (for resolved topologies) or directly from taxon sets
(for star quartets):

.. code-block:: python

   from phylozoo.core.quartet import Quartet
   from phylozoo.core.split import Split

   # Create resolved quartets using splits
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))  # AB|CD topology
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))  # AC|BD topology
   q3 = Quartet(Split({"A", "D"}, {"B", "C"}))  # AD|BC topology

   # Create star quartet (unresolved)
   star = Quartet({"A", "B", "C", "D"})

The constructor automatically determines whether the quartet is resolved or star
based on the input type. For resolved quartets, the split must be a 2|2 split (non-trivial)
on exactly 4 taxa.

Accessing Quartet Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Quartets provide access to their taxa and underlying split:

.. code-block:: python

   # Access taxa
   taxa = q1.taxa   # frozenset({"A", "B", "C", "D"})

   # Access underlying split (for resolved quartets only)
   split = q1.split  # Split object or None for star quartets

**Resolution status**

The :meth:`~phylozoo.core.quartet.base.Quartet.is_resolved` and :meth:`~phylozoo.core.quartet.base.Quartet.is_star` methods indicate whether the quartet is resolved or a star:

.. code-block:: python

   is_resolved = q1.is_resolved()   # True for resolved quartets
   is_star = star.is_star()         # True for star quartets

**Circular orderings**

Quartets can be queried for circular orderings that are congruent with their topology, using the :attr:`~phylozoo.core.quartet.base.Quartet.circular_orderings` property:

.. code-block:: python

   # Get circular orderings for a resolved quartet
   orderings = q1.circular_orderings  # Returns 2 orderings for resolved quartets

   # Get circular orderings for a star quartet
   star_orderings = star.circular_orderings  # Returns all 3 orderings for star quartets

For a resolved quartet with split {a,b}|{c,d}, this returns the two circular orderings where a and b are neighbors and c and d are neighbors. For a star tree, this returns all three circular orderings.

Quartet Operations
------------------

Quartets support various operations essential for phylogenetic analysis.

**Equality and Comparison**

Quartets implement equality comparison and hashing, making them suitable for use in sets and dictionaries:

.. code-block:: python

   # Equality comparison
   quartet_copy = Quartet(Split({"A", "B"}, {"C", "D"}))
   are_equal = q1 == quartet_copy  # True

   # Quartets can be used as set elements or dictionary keys
   quartet_set = {q1, q2, q3}      # Efficient deduplication

Quartets with the same topology (whether resolved with the same split or both star quartets on the same taxa) are considered equal.

**Network Conversion**

Resolved quartets can be converted to phylogenetic networks using the :meth:`~phylozoo.core.quartet.base.Quartet.to_network` method:

.. code-block:: python

   # Convert resolved quartet to network representation
   network = q1.to_network()  # Returns a SemiDirectedPhyNetwork

This operation creates a small phylogenetic network containing just the quartet topology, useful for visualization and further analysis.

**Copy**

The :meth:`~phylozoo.core.quartet.base.Quartet.copy` method returns an independent copy of the quartet:

.. code-block:: python

   quartet_copy = q1.copy()

See Also
--------

- :doc:`API Reference <../../../api/core/quartets>` - Complete function signatures and detailed examples
- :doc:`Quartet Profiles <quartet_profile>` - Probability distributions over quartet topologies
- :doc:`Quartet Profile Sets <quartet_profile_set>` - Collections of quartet profiles
- :doc:`Splits <../splits/overview>` - Split-based representations
