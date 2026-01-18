Quartets
========

Quartets are four-taxon unrooted trees and are fundamental building blocks for 
network inference :cite:`PhyloZoo2024`. The quartet module provides classes for 
working with quartets, quartet profiles, and quartet profile sets.

Quartet Basics
--------------

A ``Quartet`` represents an unrooted tree on four taxa. There are three possible 
quartet topologies:

.. code-block:: python

   from phylozoo import Quartet
   from phylozoo.core.split import Split
   
   # Create resolved quartets
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))  # AB|CD
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))  # AC|BD
   q3 = Quartet(Split({"A", "D"}, {"B", "C"}))  # AD|BC
   
   # Create star quartet (unresolved)
   star = Quartet({"A", "B", "C", "D"})
   
   # Check properties
   is_resolved = q1.is_resolved()
   is_star = star.is_star()
   are_compatible = q1.is_compatible(q2)

Quartet Profiles
----------------

A ``QuartetProfile`` represents a probability distribution over the three possible 
quartet topologies for a given four-taxon set:

.. code-block:: python

   from phylozoo import QuartetProfile
   
   # Create a profile with weights
   profile = QuartetProfile(
       quartet1=Quartet(Split({"A", "B"}, {"C", "D"})),
       weight1=0.7,
       quartet2=Quartet(Split({"A", "C"}, {"B", "D"})),
       weight2=0.3
   )
   
   # Access weights
   w1 = profile.weight1  # 0.7
   w2 = profile.weight2  # 0.3
   w3 = profile.weight3  # 0.0 (implicit)
   
   # Access quartets
   for quartet, weight in profile:
       print(f"{quartet}: {weight}")

Quartet Profile Sets
--------------------

A ``QuartetProfileSet`` is a collection of quartet profiles, each with its own 
weight:

.. code-block:: python

   from phylozoo import QuartetProfileSet
   
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

Complete Workflow Example
--------------------------

Here's a complete example showing how to work with quartets:

.. code-block:: python

   from phylozoo import Quartet, QuartetProfile, QuartetProfileSet
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
   
   # Use distances for network inference
   # (see :doc:`Inference <inference>` for details)

Available Functions
-------------------

**Classes:**

* **Quartet** - Immutable quartet datatype representing an unrooted tree on 4 taxa. 
  Can be resolved (2|2 split) or unresolved (star tree). Supports compatibility 
  checking and circular ordering extraction. See :class:`phylozoo.core.quartet.Quartet` 
  for full API.

* **QuartetProfile** - Immutable profile for quartets on the same 4-taxon set with 
  weights. Groups multiple quartet topologies with associated probabilities. Supports 
  iteration and weight access. See :class:`phylozoo.core.quartet.QuartetProfile` for 
  full API.

* **QuartetProfileSet** - Immutable collection of quartet profiles with two-level 
  weights. Each profile has a weight, and each quartet within a profile has a weight. 
  Can be created from profiles or quartets (auto-grouped). See 
  :class:`phylozoo.core.quartet.QuartetProfileSet` for full API.

**Quartet Methods:**

* **is_resolved()** - Check if quartet is resolved (has a split). Returns boolean.
* **is_star()** - Check if quartet is unresolved (star tree). Returns boolean.
* **is_compatible(other)** - Check compatibility with another quartet. Returns boolean.
* **taxa** - Frozen set of 4 taxon labels.

**QuartetProfile Methods:**

* **weight1, weight2, weight3** - Weights for the three possible quartet topologies.
* **quartet1, quartet2, quartet3** - The three possible quartet topologies (None if 
  weight is 0).
* **taxa** - Frozen set of 4 taxon labels.

**QuartetProfileSet Methods:**

* **taxa** - Frozen set of all taxa across all profiles.
* **profiles** - Read-only mapping of taxa sets to (profile, weight) tuples.
* **get_profile(taxa)** - Get profile for a specific 4-taxon set. Returns 
  QuartetProfile or None.

**Distance Functions** (``phylozoo.core.quartet.qdistance``):

* **quartet_distance(profileset, rho)** - Compute distance matrix from quartet 
  profile set using rho vector. The rho vector specifies distance contributions for 
  different quartet topologies. Returns DistanceMatrix. See 
  :func:`phylozoo.core.quartet.qdistance.quartet_distance`.

* **quartet_distance_with_partition(profileset, rho, partition)** - Compute distance 
  matrix with additional partition information. Advanced function for partitioned 
  analyses. See :func:`phylozoo.core.quartet.qdistance.quartet_distance_with_partition`.

.. note::
   Quartet profiles and profile sets are used extensively in network inference 
   algorithms. See :doc:`Inference <inference>` for how to use quartets for network 
   reconstruction.

.. tip::
   When creating QuartetProfileSet from quartets, quartets are automatically grouped 
   by their 4-taxon sets into profiles. This makes it easy to work with collections 
   of quartets.

.. seealso::
   For network inference using quartets, see :doc:`Inference <inference>`. For 
   extracting quartets from networks, see :doc:`Networks (Advanced) <networks/advanced>`.
