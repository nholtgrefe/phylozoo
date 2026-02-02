quartet
=======

The quartet module provides classes for representing and working with quartets, which
are unrooted trees on 4 taxa. A quartet can either be resolved (with a single
non-trivial split :math:`\{a,b\} | \{c,d\}`) or unresolved (a star tree).

Main Classes
------------

.. automodule:: phylozoo.core.quartet.base
   :members:
   :show-inheritance:

Quartet Profiles
----------------

A quartet profile groups multiple quartets on the same 4-taxon set, each with an
associated weight representing the relative importance or frequency of each quartet
topology.

.. automodule:: phylozoo.core.quartet.qprofile
   :members:
   :show-inheritance:

Quartet Profile Sets
--------------------

A quartet profile set is a collection of quartet profiles, typically used for network
reconstruction algorithms such as Squirrel.

.. automodule:: phylozoo.core.quartet.qprofileset
   :members:
   :show-inheritance:

Distance Computation
--------------------

Compute distance matrices from quartet profiles using quartet distance metrics.

.. automodule:: phylozoo.core.quartet.qdistance
   :members:
   :show-inheritance:
