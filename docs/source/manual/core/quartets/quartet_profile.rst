Quartet Profiles
================

Working with Quartet Profiles
------------------------------

The :mod:`phylozoo.core.quartet` module provides the :class:`~phylozoo.core.quartet.qprofile.QuartetProfile` class,
which represents a probability distribution over the three possible quartet topologies
for a specific four-taxon set. The weights in a quartet profile always sum to 1.0
(within a small tolerance): if no weights are given, each quartet gets equal weight
:math:`1/k` (where :math:`k` is the number of quartets); if weights are given, they
must sum to 1.0. Quartet profiles are essential for modeling
uncertainty or mixed phylogenetic signals in quartet-based analyses.

Creating Quartet Profiles
^^^^^^^^^^^^^^^^^^^^^^^^^

Quartet profiles can be created by specifying quartets and their corresponding weights:

.. code-block:: python

   from phylozoo.core.quartet import QuartetProfile, Quartet
   from phylozoo.core.split import Split

   # Create quartets
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))

   # Create a profile using a dictionary (weights must sum to 1.0)
   profile = QuartetProfile({q1: 0.7, q2: 0.3})

   # Create a profile from a list of (quartet, weight) tuples (weights must sum to 1.0)
   full_profile = QuartetProfile([
       (q1, 0.5),
       (q2, 0.3),
       (Quartet(Split({"A", "D"}, {"B", "C"})), 0.2)
   ])

   # Create a profile from a list of quartets (equal weight 1/k each, so total 1.0)
   equal_profile = QuartetProfile([q1, q2])  # Each gets weight 0.5

All quartets in a profile must have the same four taxa.

Accessing Profile Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Basic properties**

Quartet profiles provide access to their taxa, quartet weights, and induced split:

.. code-block:: python

   # Access taxa
   taxa = profile.taxa    # frozenset({"A", "B", "C", "D"})

   # Get mapping of all quartets to their weights
   quartets_map = profile.quartets  # Mapping[Quartet, float]

   # Get weight for a specific quartet
   weight = profile.get_weight(some_quartet)  # Returns weight or 0.0

   # Get split if profile has a single quartet
   split = profile.split  # Split object or None

**Circular orderings**

The :attr:`~phylozoo.core.quartet.qprofile.QuartetProfile.circular_orderings` attribute
returns circular orderings that are compatible with all quartets in the profile:

.. code-block:: python

   # Get circular orderings congruent with all quartets
   orderings = profile.circular_orderings  # frozenset[CircularOrdering] or None

**Resolution status**

The :meth:`~phylozoo.core.quartet.qprofile.QuartetProfile.is_resolved` method checks whether the profile is resolved, meaning all quartets in the profile are resolved:

.. code-block:: python

   # Check if profile is resolved (all quartets resolved)
   resolved = profile.is_resolved()  # True if all quartets are resolved

**Triviality**

The :meth:`~phylozoo.core.quartet.qprofile.QuartetProfile.is_trivial` method checks whether the profile contains only one quartet:

.. code-block:: python

   # Check if profile is trivial (only one quartet)
   trivial = profile.is_trivial()  # True if profile has only one quartet

See Also
--------

- :doc:`API Reference <../../../api/core/quartets>` - Complete function signatures and detailed examples
- :doc:`Quartets <quartet>` - Individual quartet topologies
- :doc:`Quartet Profile Sets <quartet_profile_set>` - Collections of quartet profiles
- :doc:`Distance Matrices <../distance>` - Computing distances from quartet profiles
