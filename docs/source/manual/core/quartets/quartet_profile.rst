Quartet Profiles
================

The :mod:`phylozoo.core.quartet` module provides the :class:`QuartetProfile` class,
which represents a probability distribution over the three possible quartet topologies
for a specific four-taxon set. Quartet profiles are essential for modeling uncertainty
or mixed phylogenetic signals in quartet-based analyses.

All classes and functions on this page can be imported from the core quartet module:

.. code-block:: python

   from phylozoo.core.quartet import *
   # or directly
   from phylozoo.core.quartet import QuartetProfile

Working with Quartet Profiles
-----------------------------

Quartet profiles extend individual quartets by representing uncertainty or mixed
support across the three possible quartet topologies. They are fundamental to
statistical phylogenetic inference methods that need to model confidence or
frequency distributions over quartet relationships.

.. note::
   :class: dropdown

   **Implementation details**

   Quartet profiles are optimized for statistical phylogenetic analysis:

   - Maintain normalized probability distributions automatically
   - Store only non-zero weight quartets to save memory
   - Provide efficient iteration over supported topologies
   - Immutable design ensures data integrity in analysis pipelines

   For implementation details, see :mod:`src/phylozoo/core/quartet/qprofile.py`.

Creating Quartet Profiles
^^^^^^^^^^^^^^^^^^^^^^^^^

Quartet profiles can be created by specifying quartets and their corresponding weights:

.. code-block:: python

   from phylozoo.core.quartet import QuartetProfile, Quartet
   from phylozoo.core.primitives import Split

   # Create quartets
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))

   # Create a profile using a dictionary
   profile = QuartetProfile({q1: 0.7, q2: 0.3})

   # Create a profile from a list of (quartet, weight) tuples
   full_profile = QuartetProfile([
       (q1, 0.5),
       (q2, 0.3),
       (Quartet(Split({"A", "D"}, {"B", "C"})), 0.2)
   ])

   # Create a profile from a list of quartets (equal weights)
   equal_profile = QuartetProfile([q1, q2])  # Each gets weight 1.0

Weights are automatically normalized and must sum to 1.0. Profiles can have 1, 2,
or 3 supported quartet topologies.

Accessing Profile Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quartet profiles provide comprehensive access to their probability distribution:

.. code-block:: python

   # Access taxa
   taxa = profile.taxa    # frozenset({"A", "B", "C", "D"})

   # Get mapping of all quartets to their weights
   quartets_map = profile.quartets  # Mapping[Quartet, float]

   # Get total weight (should be 1.0 for normalized profiles)
   total = profile.total_weight

   # Get weight for a specific quartet
   weight = profile.get_weight(some_quartet)  # Returns weight or 0.0

   # Get split if profile has single quartet
   split = profile.split  # Split object or None

   # Get circular orderings
   orderings = profile.circular_orderings  # frozenset[CircularOrdering] or None

Profile Operations
^^^^^^^^^^^^^^^^^^

Quartet profiles support iteration and other operations essential for analysis.

**Iteration Over Topologies**

Quartet profiles can be iterated over to access all supported quartet topologies and their weights:

.. code-block:: python

   # Iterate over supported quartets and their weights
   for quartet, weight in profile:
       print(f"Topology: {quartet}, Probability: {weight}")

   # Only non-zero weight topologies are included in iteration
   supported_topologies = list(profile)  # List of (quartet, weight) tuples

   # Get profile length (number of supported topologies)
   num_topologies = len(profile)

Iteration provides an efficient way to process all quartet topologies that have non-zero probability in the profile.

**Weight Queries**

You can query the weight of specific quartet topologies within a profile:

.. code-block:: python

   # Get weight of a specific quartet
   weight = profile.get_weight(some_quartet)  # Returns weight or 0.0

   # Check if quartet is supported (has non-zero weight)
   is_supported = some_quartet in profile

These operations allow you to check the probability distribution for specific quartet relationships.

**Resolution and Properties**

Quartet profiles provide methods to check their mathematical properties:

.. code-block:: python

   # Check if profile represents resolved relationships
   resolved = profile.is_resolved()  # True if all quartets are resolved

   # Check if profile is trivial (all weight on one topology)
   trivial = profile.is_trivial()

   # Get single split if profile has only one topology
   split = profile.split  # Split or None

   # Get circular orderings if applicable
   orderings = profile.circular_orderings  # frozenset[CircularOrdering] or None

**Mathematical Properties**

Quartet profiles implement equality comparison and hashing:

.. code-block:: python

   # Profiles are immutable and hashable
   profile_copy = QuartetProfile(profile.quartets)  # Copy from quartets mapping
   are_equal = profile == profile_copy  # True

   # Can be used as dictionary keys or set elements
   profile_dict = {profile: "data"}
   profile_set = {profile, profile_copy}  # Contains only one element

   # String representations
   profile_str = str(profile)
   profile_repr = repr(profile)

Profiles with identical quartet topologies and weights are considered equal, enabling efficient storage and comparison.

See Also
--------

- :doc:`API Reference <../../../api/core/quartet>` - Complete function signatures and detailed examples
- :doc:`Quartets <quartet>` - Individual quartet topologies
- :doc:`Quartet Profile Sets <quartet_profile_set>` - Collections of quartet profiles
- :doc:`Distance Matrices <../distance>` - Computing distances from quartet profiles
