Triplets
========

Working with Triplets
---------------------

The :mod:`phylozoo.core.triplet` module provides the :class:`~phylozoo.core.triplet.base.Triplet` class, which
represents a rooted tree topology on three taxa. Triplets are the fundamental units
of rooted phylogenetic relationships in triplet-based methods, with three possible
resolved topologies (one per choice of outgroup) or an unresolved star triplet.

Creating Triplets
^^^^^^^^^^^^^^^^^

Triplets can be created from splits (for resolved topologies) or directly from taxon sets
(for star triplets):

.. code-block:: python

   from phylozoo.core.triplet import Triplet
   from phylozoo.core.split import Split

   # Create resolved triplets using trivial splits a|bc
   t1 = Triplet(Split({"A"}, {"B", "C"}))  # A|BC topology (A is outgroup)
   t2 = Triplet(Split({"B"}, {"A", "C"}))  # B|AC topology (B is outgroup)
   t3 = Triplet(Split({"C"}, {"A", "B"}))  # C|AB topology (C is outgroup)

   # Create star triplet (unresolved)
   star = Triplet({"A", "B", "C"})

The constructor automatically determines whether the triplet is resolved or star
based on the input type. For resolved triplets, the split must be a trivial 1|2 split
on exactly 3 taxa.

Accessing Triplet Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Triplets provide access to their taxa and underlying split:

.. code-block:: python

   # Access taxa
   taxa = t1.taxa   # frozenset({"A", "B", "C"})

   # Access underlying split (for resolved triplets only)
   split = t1.split  # Split object or None for star triplets

**Resolution status**

The :meth:`~phylozoo.core.triplet.base.Triplet.is_resolved` and :meth:`~phylozoo.core.triplet.base.Triplet.is_star` methods indicate whether the triplet is resolved or a star:

.. code-block:: python

   is_resolved = t1.is_resolved()   # True for resolved triplets
   is_star = star.is_star()         # True for star triplets

**Outgroup and cherry**

Resolved triplets expose their outgroup and cherry through the :attr:`~phylozoo.core.triplet.base.Triplet.outgroup` and :attr:`~phylozoo.core.triplet.base.Triplet.cherry` properties:

.. code-block:: python

   # For a resolved triplet a|bc:
   outgroup = t1.outgroup   # frozenset({"A"})   - the lone child of the root
   cherry = t1.cherry       # frozenset({"B", "C"}) - the two siblings under an internal node

   # For a star triplet:
   star.outgroup is None    # True
   star.cherry is None      # True

Triplet Operations
------------------

Triplets support various operations essential for phylogenetic analysis.

**Equality and Comparison**

Triplets implement equality comparison and hashing, making them suitable for use in sets and dictionaries:

.. code-block:: python

   # Equality comparison
   triplet_copy = Triplet(Split({"A"}, {"B", "C"}))
   are_equal = t1 == triplet_copy  # True

   # Triplets can be used as set elements or dictionary keys
   triplet_set = {t1, t2, t3}      # Efficient deduplication

Triplets with the same topology (whether resolved with the same split or both star triplets on the same taxa) are considered equal.

**Network Conversion**

Triplets can be converted to phylogenetic networks using the :meth:`~phylozoo.core.triplet.base.Triplet.to_network` method:

.. code-block:: python

   # Convert triplet to rooted network representation
   network = t1.to_network()  # Returns a DirectedPhyNetwork

This operation creates a small rooted phylogenetic network containing just the triplet topology, useful for visualization and further analysis.

**Copy**

The :meth:`~phylozoo.core.triplet.base.Triplet.copy` method returns an independent copy of the triplet:

.. code-block:: python

   triplet_copy = t1.copy()

See Also
--------

- :doc:`API Reference <../../../api/core/triplets>` - Complete function signatures and detailed examples
- :doc:`Triplet Profiles <triplet_profile>` - Probability distributions over triplet topologies
- :doc:`Triplet Profile Sets <triplet_profile_set>` - Collections of triplet profiles
- :doc:`Splits <../splits/overview>` - Split-based representations
- :doc:`Quartets <../quartets/quartet>` - Unrooted four-taxon analogues
