Quartet Profile Sets
=====================

The :mod:`phylozoo.core.quartet` module provides the :class:`QuartetProfileSet` class,
which represents a weighted collection of quartet profiles covering multiple four-taxon
subsets. This is the primary data structure for quartet-based phylogenetic analysis
and network inference algorithms.

All classes and functions on this page can be imported from the core quartet module:

.. code-block:: python

   from phylozoo.core.quartet import *
   # or directly
   from phylozoo.core.quartet import QuartetProfileSet

Working with Quartet Profile Sets
----------------------------------

Quartet profile sets are the fundamental data structure for quartet-based phylogenetic
methods. They manage collections of quartet profiles, automatically organizing them
by four-taxon sets and providing efficient access for analysis and inference algorithms.

.. note::
   :class: dropdown

   **Implementation details**

   Quartet profile sets are optimized for large-scale phylogenetic analysis:

   - Efficient indexing by four-taxon sets for fast profile lookup
   - Automatic quartet aggregation when created from individual quartets
   - Memory-efficient storage of sparse quartet systems
   - Immutable design ensuring data integrity in analysis pipelines

   For implementation details, see :mod:`src/phylozoo/core/quartet/qprofileset.py`.

Creating Quartet Profile Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quartet profile sets can be created from existing profiles or directly from quartets:

**From Quartet Profiles**

.. code-block:: python

   from phylozoo.core.quartet import QuartetProfileSet, QuartetProfile, Quartet
   from phylozoo.core.primitives import Split

   # Create individual profiles
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))

   profile1 = QuartetProfile({q1: 0.8})
   profile2 = QuartetProfile({q2: 0.6})

   # Create profile set with profile-weight tuples
   profile_set = QuartetProfileSet([
       (profile1, 0.5),
       (profile2, 0.5)
   ])

**From Individual Quartets**

.. code-block:: python

   # Create from quartets (automatically grouped by taxa)
   quartets = [
       Quartet(Split({"A", "B"}, {"C", "D"})),
       Quartet(Split({"A", "C"}, {"B", "D"})),
       Quartet(Split({"A", "B"}, {"C", "E"})),  # Different 4-taxon set
   ]
   profile_set = QuartetProfileSet(quartets)

When created from quartets, they are automatically grouped by their four-taxon sets
and converted to profiles.

Accessing Profile Set Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quartet profile sets provide comprehensive access to their structure and contents:

.. code-block:: python

   # Basic properties
   total_taxa = profile_set.taxa        # frozenset of all taxa
   num_profiles = len(profile_set)      # Number of profiles

   # Check density and resolution properties
   is_dense = profile_set.is_dense      # True if has all possible 4-taxon combinations
   all_resolved = profile_set.is_all_resolved  # True if all profiles are resolved
   max_len = profile_set.max_profile_len       # Maximum quartets per profile

   # Access individual profiles
   profile = profile_set.get_profile(frozenset({"A", "B", "C", "D"}))
   # Returns QuartetProfile or None

   profile_weight = profile_set.get_profile_weight(frozenset({"A", "B", "C", "D"}))
   # Returns float or None

   has_profile = profile_set.has_profile(frozenset({"A", "B", "C", "D"}))

   # Access all profiles (read-only mapping)
   all_profiles = profile_set.profiles  # Dict[frozenset, (QuartetProfile, float)]

Profile Set Operations
^^^^^^^^^^^^^^^^^^^^^^

Quartet profile sets support various operations for analysis and processing.

**Iteration and Membership**

Quartet profile sets can be iterated over and queried for membership:

.. code-block:: python

   # Iterate over all profiles in the set
   for profile, weight in profile_set:
       print(f"Profile weight: {weight}")

   # Check if a specific 4-taxon set has a profile
   has_four_taxa = frozenset({"A", "B", "C", "D"}) in profile_set

   # Get all 4-taxon sets that have profiles
   all_taxon_sets = list(profile_set.all_profile_taxon_sets())

These operations allow you to traverse and query the quartet profile collection efficiently.

**Distance Matrix Computation**

The primary operation on quartet profile sets is computing evolutionary distance matrices:

.. code-block:: python

   from phylozoo.core.quartet.qdistance import quartet_distance

   # Compute distance matrix using rho vector
   # rho = (rho_c, rho_s, rho_a, rho_o) - quartet topology weights
   rho = (1.0, 0.5, 0.8, 0.3)
   distance_matrix = quartet_distance(profile_set, rho)

The rho vector specifies how different quartet topologies contribute to pairwise distances between taxa. This enables flexible distance computations for various phylogenetic applications. The quartet distance computation is fundamental for quartet-based phylogenetic inference methods.

**Advanced Distance Computation**

For specialized analyses, quartet profile sets support distance computation with additional partitioning:

.. code-block:: python

   from phylozoo.core.quartet.qdistance import quartet_distance_with_partition
   from phylozoo.core.primitives import Partition

   # Compute distance matrix with partition information
   partition = Partition([{"A", "B"}, {"C", "D", "E"}])
   distance_matrix = quartet_distance_with_partition(profile_set, rho, partition)

This advanced method incorporates additional partitioning information for specialized analyses, allowing you to account for known groupings or constraints in the distance computation.

See Also
--------

- :doc:`API Reference <../../../api/core/quartet>` - Complete function signatures and detailed examples
- :doc:`Quartets <quartet>` - Individual quartet topologies
- :doc:`Quartet Profiles <quartet_profile>` - Probability distributions over quartet topologies
- :doc:`Distance Matrices <../distance>` - Distance matrix computations