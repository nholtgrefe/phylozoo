Quartets
========

The :mod:`phylozoo.core.quartet` module provides the :class:`Quartet` class, which
represents an unrooted tree topology on four taxa. Quartets are the fundamental units
of phylogenetic relationships in quartet-based methods, with three possible resolved
topologies or an unresolved star quartet.

All classes and functions on this page can be imported from the core quartet module:

.. code-block:: python

   from phylozoo.core.quartet import *
   # or directly
   from phylozoo.core.quartet import Quartet

Working with Quartets
----------------------

Quartets represent the smallest phylogenetic relationships that can exist between
four taxa. Each quartet encodes either a resolved tree topology (one of three possible
bipartitions) or an unresolved star relationship where all four taxa are equally
related.

.. note::
   :class: dropdown

   **Implementation details**

   Quartet representation is optimized for phylogenetic operations:

   - Resolved quartets store their topology as a split for efficient compatibility checking
   - Star quartets use a frozenset representation for the four taxa
   - Equality and hashing are based on canonical representations
   - Immutable design ensures data integrity in analysis pipelines

   For implementation details, see :mod:`src/phylozoo/core/quartet/base.py`.

Creating Quartets
^^^^^^^^^^^^^^^^^

Quartets can be created from splits (for resolved topologies) or directly from taxon sets
(for star quartets):

.. code-block:: python

   from phylozoo.core.quartet import Quartet
   from phylozoo.core.primitives import Split

   # Create resolved quartets using splits
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))  # AB|CD topology
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))  # AC|BD topology
   q3 = Quartet(Split({"A", "D"}, {"B", "C"}))  # AD|BC topology

   # Create star quartet (unresolved)
   star = Quartet({"A", "B", "C", "D"})

The constructor automatically determines whether the quartet is resolved or star
based on the input type.

Accessing Quartet Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quartets provide methods to query their structure and relationships:

.. code-block:: python

   # Check resolution status
   is_resolved = q1.is_resolved()    # True for resolved quartets
   is_star = star.is_star()          # True for star quartets

   # Access taxa
   taxa = q1.taxa                    # frozenset({"A", "B", "C", "D"})

   # Access underlying split (for resolved quartets)
   split = q1.split                  # Split object

   # Get circular orderings (for resolved quartets)
   orderings = q1.circular_orderings  # frozenset of CircularOrdering objects

   # Create a copy
   quartet_copy = q1.copy()

Quartet Operations
^^^^^^^^^^^^^^^^^^

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

**String Representations**

Quartets provide both human-readable and programmatic string representations:

.. code-block:: python

   # Human-readable string representation
   quartet_str = str(q1)           # Shows topology in readable format

   # Python representation (for debugging)
   quartet_repr = repr(q1)         # Shows constructor call

These representations help with debugging and logging quartet-based analyses.

**Network Conversion**

Resolved quartets can be converted to phylogenetic networks:

.. code-block:: python

   # Convert resolved quartet to network representation
   network = q1.to_network()       # Returns SemiDirectedPhyNetwork

This operation creates a small phylogenetic network containing just the quartet topology, useful for visualization and further analysis. Star quartets cannot be converted to networks as they represent unresolved relationships.

See Also
--------

- :doc:`API Reference <../../../api/core/quartet>` - Complete function signatures and detailed examples
- :doc:`Quartet Profiles <quartet_profile>` - Probability distributions over quartet topologies
- :doc:`Quartet Profile Sets <quartet_profile_set>` - Collections of quartet profiles
- :doc:`Splits <../splits/overview>` - Split-based representations
