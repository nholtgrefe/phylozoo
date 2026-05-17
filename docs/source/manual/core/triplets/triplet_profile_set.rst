Triplet Profile Sets
====================

Working with Triplet Profile Sets
---------------------------------

The :mod:`phylozoo.core.triplet` module provides the :class:`~phylozoo.core.triplet.tprofileset.TripletProfileSet` class,
which represents a weighted collection of triplet profiles covering multiple three-taxon
sets.

Note that this allows for two-level weights: the weight of a profile and the weight of a triplet within a profile.

Creating Triplet Profile Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Triplet profile sets can be created from existing profiles or directly from triplets:

**From Triplet Profiles**

.. code-block:: python

   from phylozoo.core.triplet import TripletProfileSet, TripletProfile, Triplet
   from phylozoo.core.split import Split

   # Create individual profiles
   t1 = Triplet(Split({"A"}, {"B", "C"}))
   t2 = Triplet(Split({"B"}, {"A", "C"}))

   profile1 = TripletProfile({t1: 1.0})  # Single triplet: weight must be 1.0
   profile2 = TripletProfile({t2: 1.0})

   # Create profile set with profile-weight tuples
   profile_set = TripletProfileSet([
       (profile1, 0.5),
       (profile2, 0.5)
   ])

**From Individual Triplets**

.. code-block:: python

   # Create from triplets (automatically grouped by taxa)
   triplets = [
       Triplet(Split({"A"}, {"B", "C"})),
       Triplet(Split({"B"}, {"A", "C"})),
       Triplet(Split({"A"}, {"B", "D"})),  # Different 3-taxon set
   ]
   profile_set = TripletProfileSet(profiles=triplets)

When created from triplets, they are automatically grouped by their three-taxon sets
and converted to profiles. For each 3-taxon set, all triplets on that set form a
single :class:`~phylozoo.core.triplet.tprofile.TripletProfile` in which every triplet
receives equal weight :math:`1/k` (where :math:`k` is the number of triplets for that
taxa set); the resulting profile in the set has default profile weight 1.0.

If you need non-uniform triplet weights within a profile, construct a
:class:`~phylozoo.core.triplet.tprofile.TripletProfile` explicitly (using a dictionary
or list of ``(Triplet, weight)`` pairs that sum to 1.0) and pass that
``TripletProfile`` to :class:`~phylozoo.core.triplet.tprofileset.TripletProfileSet`,
optionally together with a separate profile weight.

**Specifying Total Taxa**

You can also specify the total set of taxa, which allows including taxa that don't
appear in any profile:

.. code-block:: python

   # Create profile set with explicit taxa set
   profile_set = TripletProfileSet(
       profiles=[profile1, profile2],
       taxa=frozenset({"A", "B", "C", "D", "E"})
   )

Accessing Profile Set Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Triplet profile sets provide comprehensive access to their structure and contents:

.. code-block:: python

   # Basic properties
   total_taxa = profile_set.taxa        # frozenset of all taxa
   num_profiles = len(profile_set)      # Number of profiles

   # Check maximum triplets per profile
   max_len = profile_set.max_profile_len

   # Access individual profiles
   profile = profile_set.get_profile(frozenset({"A", "B", "C"}))
   # Returns TripletProfile or None

   profile_weight = profile_set.get_profile_weight(frozenset({"A", "B", "C"}))
   # Returns float or None

   has_profile = profile_set.has_profile(frozenset({"A", "B", "C"}))

   # Access all profiles (read-only mapping)
   all_profiles = profile_set.profiles  # Dict[frozenset, (TripletProfile, float)]


**Density**

The :attr:`~phylozoo.core.triplet.tprofileset.TripletProfileSet.is_dense` property checks if the triplet profile set is dense, meaning it has a profile for every possible 3-taxon combination.

.. code-block:: python

   is_dense = profile_set.is_dense  # True if has all possible 3-taxon combinations


**Resolution status**

The :attr:`~phylozoo.core.triplet.tprofileset.TripletProfileSet.is_all_resolved` property checks if all profiles in the set are resolved, meaning all triplets in the profile are resolved.

.. code-block:: python

   is_all_resolved = profile_set.is_all_resolved  # True if all profiles are resolved


See Also
--------

- :doc:`API Reference <../../../api/core/triplets>` - Complete function signatures and detailed examples
- :doc:`Triplets <triplet>` - Individual rooted triplet topologies
- :doc:`Triplet Profiles <triplet_profile>` - Sets of triplets on the same 3-taxon set with weights
- :doc:`Quartet Profile Sets <../quartets/quartet_profile_set>` - Unrooted four-taxon analogues
