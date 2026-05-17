Overview
========

The :mod:`phylozoo.core.triplet` module provides data structures for working with
triplets, which are three-taxon rooted trees that serve as building blocks for
many concepts in phylogenetic analysis.

All classes and functions on this page can be imported from the core triplet module:

.. code-block:: python

   from phylozoo.core.triplet import *
   # or directly
   from phylozoo.core.triplet import Triplet, TripletProfile, TripletProfileSet

Classes
-------

- :doc:`Triplet <triplet>` - Single rooted triplet topologies
- :doc:`Triplet Profile <triplet_profile>` - Set of triplets on the same 3-taxon set with weights
- :doc:`Triplet Profile Set <triplet_profile_set>` - Collections of triplet profiles

See Also
--------

- :doc:`API Reference <../../../api/core/triplets>` - Complete function signatures and detailed examples
- :doc:`Quartets <../quartets/overview>` - Unrooted four-taxon analogues
