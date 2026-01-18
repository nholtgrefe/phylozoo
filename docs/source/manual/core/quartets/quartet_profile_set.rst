QuartetProfileSet
=================

A ``QuartetProfileSet`` is a collection of quartet profiles, each with its own 
weight. It can be created from profiles or directly from quartets (which are 
automatically grouped by their 4-taxon sets).

Creating Profile Sets
---------------------

.. code-block:: python

   from phylozoo import Quartet, QuartetProfile, QuartetProfileSet
   from phylozoo.core.split import Split
   
   # Create profile set from profiles
   profiles = [
       QuartetProfile(
           Quartet(Split({"A", "B"}, {"C", "D"})), 
           weight1=0.8
       ),
       QuartetProfile(
           Quartet(Split({"A", "C"}, {"B", "D"})), 
           weight1=0.6
       ),
   ]
   
   profile_set = QuartetProfileSet(profiles, weights=[0.5, 0.5])
   
   # Or create from quartets directly (auto-grouped by taxa)
   quartets = [
       Quartet(Split({"A", "B"}, {"C", "D"})),
       Quartet(Split({"A", "C"}, {"B", "D"})),
   ]
   profile_set = QuartetProfileSet(quartets)

Working with Profile Sets
-------------------------

.. code-block:: python

   # Access properties
   taxa = profile_set.taxa  # Frozen set of all taxa
   num_profiles = len(profile_set)
   
   # Get profile for a specific 4-taxon set
   profile = profile_set.get_profile(frozenset({"A", "B", "C", "D"}))
   
   # Access profiles mapping
   profiles = profile_set.profiles  # Read-only mapping

Distance Calculations
---------------------

Quartet profiles can be used to compute distance matrices:

.. code-block:: python

   from phylozoo.core.quartet.qdistance import quartet_distance
   
   # Compute distance matrix from quartet profile set
   # rho vector: (rho_c, rho_s, rho_a, rho_o)
   rho = (1.0, 0.5, 0.8, 0.3)
   distance_matrix = quartet_distance(profile_set, rho)

The rho vector specifies how different quartet topologies contribute to distances. 
See :func:`phylozoo.core.quartet.qdistance.quartet_distance` for details.

Complete Example
-----------------

.. code-block:: python

   from phylozoo import Quartet, QuartetProfileSet
   from phylozoo.core.split import Split
   from phylozoo.core.quartet.qdistance import quartet_distance
   
   # Create quartets for multiple 4-taxon sets
   quartets = [
       Quartet(Split({"A", "B"}, {"C", "D"})),
       Quartet(Split({"A", "B"}, {"C", "E"})),
       Quartet(Split({"A", "C"}, {"D", "E"})),
   ]
   
   # Create profile set (auto-groups by taxa)
   profile_set = QuartetProfileSet(quartets)
   
   # Check properties
   print(f"Number of profiles: {len(profile_set)}")
   print(f"Total taxa: {profile_set.taxa}")
   
   # Compute distance matrix
   rho = (1.0, 0.5, 0.8, 0.3)
   distances = quartet_distance(profile_set, rho)

API Reference
-------------

**Class**: :class:`phylozoo.core.quartet.QuartetProfileSet`

**Properties:**

* **taxa** - Frozen set of all taxa across all profiles.
* **profiles** - Read-only mapping of taxa sets to (profile, weight) tuples.

**Methods:**

* **get_profile(taxa)** - Get profile for a specific 4-taxon set. Returns 
  QuartetProfile or None.

**Distance Functions** (``phylozoo.core.quartet.qdistance``):

* **quartet_distance(profileset, rho)** - Compute distance matrix from quartet 
  profile set using rho vector. Returns DistanceMatrix. See 
  :func:`phylozoo.core.quartet.qdistance.quartet_distance`.

* **quartet_distance_with_partition(profileset, rho, partition)** - Compute distance 
  matrix with additional partition information. Advanced function for partitioned 
  analyses. See :func:`phylozoo.core.quartet.qdistance.quartet_distance_with_partition`.

.. tip::
   When creating QuartetProfileSet from quartets, quartets are automatically grouped 
   by their 4-taxon sets into profiles. This makes it easy to work with collections 
   of quartets.

.. seealso::
   For single quartets, see :doc:`Quartet <quartet>`. 
   For quartet profiles, see :doc:`QuartetProfile <quartet_profile>`. 
   For network inference using quartets, see :doc:`Inference <../../inference/overview>`.
