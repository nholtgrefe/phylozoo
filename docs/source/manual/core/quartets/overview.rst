Overview
========

The :mod:`phylozoo.core.quartet` module provides data structures for working with
quartets, which are four-taxon unrooted trees that serve as building blocks for
phylogenetic network inference algorithms.

All classes and functions on this page can be imported from the core quartet module:

.. code-block:: python

   from phylozoo.core.quartet import *
   # or directly
   from phylozoo.core.quartet import Quartet, QuartetProfile, QuartetProfileSet

Classes
-------

- :doc:`Quartet <quartet>` - Single quartet topologies
- :doc:`Quartet Profile <quartet_profile>` - Probability distributions over quartet topologies
- :doc:`Quartet Profile Set <quartet_profile_set>` - Collections of quartet profiles

See Also
--------

- :doc:`API Reference <../../../api/core/quartet>` - Complete function signatures and detailed examples
- :doc:`Distance Matrices <../distance>` - Computing distances from quartet profiles
- :doc:`Inference <../../inference/overview>` - Network inference using quartets
