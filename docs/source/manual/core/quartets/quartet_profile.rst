QuartetProfile
==============

A ``QuartetProfile`` represents a probability distribution over the three possible 
quartet topologies for a given four-taxon set.

Creating Quartet Profiles
--------------------------

.. code-block:: python

   from phylozoo import Quartet, QuartetProfile
   from phylozoo.core.split import Split
   
   # Create a profile with weights
   profile = QuartetProfile(
       quartet1=Quartet(Split({"A", "B"}, {"C", "D"})),
       weight1=0.7,
       quartet2=Quartet(Split({"A", "C"}, {"B", "D"})),
       weight2=0.3
   )

Accessing Weights and Quartets
-------------------------------

.. code-block:: python

   # Access weights
   w1 = profile.weight1  # 0.7
   w2 = profile.weight2  # 0.3
   w3 = profile.weight3  # 0.0 (implicit)
   
   # Access quartets
   q1 = profile.quartet1
   q2 = profile.quartet2
   q3 = profile.quartet3  # None if weight is 0
   
   # Iterate over quartets and weights
   for quartet, weight in profile:
       print(f"{quartet}: {weight}")

API Reference
-------------

**Class**: :class:`phylozoo.core.quartet.QuartetProfile`

**Properties:**

* **weight1, weight2, weight3** - Weights for the three possible quartet topologies.
* **quartet1, quartet2, quartet3** - The three possible quartet topologies (None if 
  weight is 0).
* **taxa** - Frozen set of 4 taxon labels.

**Iteration:**

* Iterating over a profile yields ``(quartet, weight)`` tuples for non-zero weights.

.. seealso::
   For single quartets, see :doc:`Quartet <quartet>`. 
   For collections of profiles, see :doc:`QuartetProfileSet <quartet_profile_set>`.
