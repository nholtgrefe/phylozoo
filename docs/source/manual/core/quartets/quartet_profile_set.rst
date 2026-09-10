Quartet Profile Sets
====================

Working with Quartet Profile Sets
----------------------------------

The :mod:`phylozoo.core.quartet` module provides the :class:`~phylozoo.core.quartet.qprofileset.QuartetProfileSet` class,
which represents a weighted collection of quartet profiles covering multiple four-taxon
.

Note that this allows for two-level weights: the weight of a profile and the weight of a quartet within a profile.

Creating Quartet Profile Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quartet profile sets can be created from existing profiles or directly from quartets:

**From Quartet Profiles**

.. code-block:: python

   from phylozoo.core.quartet import QuartetProfileSet, QuartetProfile, Quartet
   from phylozoo.core.split import Split

   # Create individual profiles
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))

   profile1 = QuartetProfile({q1: 1.0})  # Single quartet: weight must be 1.0
   profile2 = QuartetProfile({q2: 1.0})

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
   profile_set = QuartetProfileSet(profiles=quartets)

When created from quartets, they are automatically grouped by their four-taxon sets
and converted to profiles. For each 4-taxon set, all quartets on that set form a
single :class:`~phylozoo.core.quartet.qprofile.QuartetProfile` in which every quartet
receives equal weight :math:`1/k` (where :math:`k` is the number of quartets for that
taxa set); the resulting profile in the set has default profile weight 1.0.

If you need non-uniform quartet weights within a profile, construct a
:class:`~phylozoo.core.quartet.qprofile.QuartetProfile` explicitly (using a dictionary
or list of ``(Quartet, weight)`` pairs that sum to 1.0) and pass that
``QuartetProfile`` to :class:`~phylozoo.core.quartet.qprofileset.QuartetProfileSet`,
optionally together with a separate profile weight.

**Specifying Total Taxa**

You can also specify the total set of taxa, which allows including taxa that don't
appear in any profile:

.. code-block:: python

   # Create profile set with explicit taxa set
   profile_set = QuartetProfileSet(
       profiles=[profile1, profile2],
       taxa=frozenset({"A", "B", "C", "D", "E", "F"})
   )

Accessing Profile Set Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Quartet profile sets provide comprehensive access to their structure and contents:

.. code-block:: python

   # Basic properties
   total_taxa = profile_set.taxa        # frozenset of all taxa
   num_profiles = len(profile_set)      # Number of profiles

   # Check maximum quartets per profile
   max_len = profile_set.max_profile_len

   # Access individual profiles
   profile = profile_set.get_profile(frozenset({"A", "B", "C", "D"}))
   # Returns QuartetProfile or None

   profile_weight = profile_set.get_profile_weight(frozenset({"A", "B", "C", "D"}))
   # Returns float or None

   has_profile = profile_set.has_profile(frozenset({"A", "B", "C", "D"}))

   # Access all profiles (read-only mapping)
   all_profiles = profile_set.profiles  # Dict[frozenset, (QuartetProfile, float)]


**Density**

The :attr:`~phylozoo.core.quartet.qprofileset.QuartetProfileSet.is_dense` property checks if the quartet profile set is dense, meaning it has a profile for every possible 4-taxon combination.

.. code-block:: python

   is_dense = profile_set.is_dense  # True if has all possible 4-taxon combinations


**Resolution status**

The :attr:`~phylozoo.core.quartet.qprofileset.QuartetProfileSet.is_all_resolved` property checks if all profiles in the set are resolved, meaning all quartets in the profile are resolved.

.. code-block:: python

   is_all_resolved = profile_set.is_all_resolved  # True if all profiles are resolved


Quartet Distance Computation
----------------------------

The quartet module provides functions for computing distance matrices from quartet profile sets.
This quartet distance metric was first defined for trees and their quartets :cite:`Rhodes2019`,
then extended to networks for the NANUQ algorithm :cite:`Allman2019`, and later further explored in various forms in :cite:`Holtgrefe2025Squirrel` (for the Squirrel algorithm), :cite:`Allman2025` (for the NANUQ+ algorithm), and :cite:`Holtgrefe2025b` (for level-2 networks).

Quartet Distance
^^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.quartet.qdistance.quartet_distance` function computes
a distance matrix from a quartet profile set using a rho vector.
The distance between two taxa is computed by aggregating contributions from all
quartet profiles, where the contribution depends on the quartet topology and the
rho vector values. The current implementation allows only dense quartet profile sets
with exactly 1 or 2 resolved quartets per 4-taxon set.

.. code-block:: python

   from phylozoo.core.quartet.qdistance import quartet_distance

   # Compute distance matrix using rho vector
   rho = (0.5, 1.0, 0.5, 1.0)  # Squirrel rho vector
   distance_matrix = quartet_distance(profile_set, rho)

The distance formula computes pairwise distances between taxa :math:`i` and :math:`j` as:

.. math::

   D(i,j) = \sum_{\substack{S \subseteq X \\ |S| = 4 \\ i,j \in S}} 2 \cdot \rho_{\text{dist}}(Q_S, i, j, \rho) + 2n - 4

where :math:`n` is the number of taxa, :math:`X` is the set of all taxa, :math:`Q_S` is the quartet profile
for the 4-taxon set :math:`S`, and :math:`\rho_{\text{dist}}` is the rho-distance function.

The rho-distance function :math:`\rho_{\text{dist}}(Q_S, i, j, \rho)` depends on the quartet profile type,
the leaves :math:`i` and :math:`j`, and the rho vector :math:`\rho = (\rho_c, \rho_s, \rho_a, \rho_o)`

*Profiles with 1 quartet (split):*

For a profile containing a single quartet with split :math:`\{a,b\} | \{c,d\}`, the rho-distance is:

.. math::

   \rho_{\text{dist}}(Q_S, i, j, \rho) = \begin{cases}
   \rho_c & \text{if } i \text{ and } j \text{ are on the same side of the split} \\
   \rho_s & \text{if } i \text{ and } j \text{ are on different sides of the split}
   \end{cases}

Note: The naming convention uses :math:`\rho_c` for "cherry" (same side) and :math:`\rho_s` for "split" (different sides).

*Profiles with 2 quartets (four-cycle):*

For a profile containing two quartets, which therefore induce a single circular ordering, the rho-distance is:

.. math::

   \rho_{\text{dist}}(Q_S, i, j, \rho) = \begin{cases}
   \rho_a & \text{if } i \text{ and } j \text{ are adjacent in the circular ordering} \\
   \rho_o & \text{if } i \text{ and } j \text{ are opposite in the circular ordering}
   \end{cases}

*Note:*
The rho vector must satisfy :math:`\rho_a \leq \rho_o` and :math:`\rho_c \leq \rho_s`.

Common rho vector values:

- **NANUQ**: :math:`(0.0, 1.0, 0.5, 1.0)` :cite:`Allman2019`, :cite:`Holtgrefe2025b`
- **Squirrel/MONAD**: :math:`(0.5, 1.0, 0.5, 1.0)` :cite:`Holtgrefe2025Squirrel`, :cite:`Allman2025`

See Also
--------

- :doc:`API Reference <../../../api/core/quartets>` - Complete function signatures and detailed examples
- :doc:`Quartets <quartet>` - Individual quartet topologies
- :doc:`Quartet Profiles <quartet_profile>` - Sets of quartets on the same 4-taxon set with weights
- :doc:`Distance Matrices <../distance>` - Distance matrix computations
