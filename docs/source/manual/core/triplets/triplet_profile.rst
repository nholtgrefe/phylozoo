Triplet Profiles
================

Working with Triplet Profiles
-----------------------------

The :mod:`phylozoo.core.triplet` module provides the :class:`~phylozoo.core.triplet.tprofile.TripletProfile` class,
which represents a probability distribution over the possible triplet topologies
for a specific three-taxon set. The weights in a triplet profile always sum to 1.0
(within a small tolerance): if no weights are given, each triplet gets equal weight
:math:`1/k` (where :math:`k` is the number of triplets); if weights are given, they
must sum to 1.0. Triplet profiles are essential for modeling
uncertainty or mixed phylogenetic signals in triplet-based analyses.

Creating Triplet Profiles
^^^^^^^^^^^^^^^^^^^^^^^^^

Triplet profiles can be created by specifying triplets and their corresponding weights:

.. code-block:: python

   from phylozoo.core.triplet import TripletProfile, Triplet
   from phylozoo.core.split import Split

   # Create triplets
   t1 = Triplet(Split({"A"}, {"B", "C"}))
   t2 = Triplet(Split({"B"}, {"A", "C"}))

   # Create a profile using a dictionary (weights must sum to 1.0)
   profile = TripletProfile({t1: 0.7, t2: 0.3})

   # Create a profile from a list of (triplet, weight) tuples (weights must sum to 1.0)
   full_profile = TripletProfile([
       (t1, 0.5),
       (t2, 0.3),
       (Triplet(Split({"C"}, {"A", "B"})), 0.2)
   ])

   # Create a profile from a list of triplets (equal weight 1/k each, so total 1.0)
   equal_profile = TripletProfile([t1, t2])  # Each gets weight 0.5

All triplets in a profile must have the same three taxa.

Accessing Profile Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Triplet profiles provide access to their taxa, triplet weights, and induced split:

.. code-block:: python

   # Access taxa
   taxa = profile.taxa    # frozenset({"A", "B", "C"})

   # Get mapping of all triplets to their weights
   triplets_map = profile.triplets  # Mapping[Triplet, float]

   # Get weight for a specific triplet
   weight = profile.get_weight(some_triplet)  # Returns weight or 0.0

   # Get split if profile has a single triplet
   split = profile.split  # Split object or None

**Resolution status**

The :meth:`~phylozoo.core.triplet.tprofile.TripletProfile.is_resolved` method checks whether the profile is resolved, meaning all triplets in the profile are resolved:

.. code-block:: python

   # Check if profile is resolved (all triplets resolved)
   resolved = profile.is_resolved()  # True if all triplets are resolved

**Triviality**

The :meth:`~phylozoo.core.triplet.tprofile.TripletProfile.is_trivial` method checks whether the profile contains only one triplet:

.. code-block:: python

   # Check if profile is trivial (only one triplet)
   trivial = profile.is_trivial()  # True if profile has only one triplet

See Also
--------

- :doc:`API Reference <../../../api/core/triplets>` - Complete function signatures and detailed examples
- :doc:`Triplets <triplet>` - Individual rooted triplet topologies
- :doc:`Triplet Profile Sets <triplet_profile_set>` - Collections of triplet profiles
- :doc:`Quartet Profiles <../quartets/quartet_profile>` - Unrooted four-taxon analogues
