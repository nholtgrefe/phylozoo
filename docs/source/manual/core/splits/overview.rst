Splits
======

The :mod:`phylozoo.core.split` module provides data structures for working with
splits, which represent bipartitions of taxa sets and form the mathematical foundation
for phylogenetic tree and network representations.

All classes and functions on this page can be imported from the core split module:

.. code-block:: python

   from phylozoo.core.split import *
   # or directly
   from phylozoo.core.split import Split, SplitSystem

Classes
-------

- :doc:`Split <split>` - Individual bipartitions
- :doc:`Split System <split_system>` - Collections of splits for phylogenetic analysis

See Also
--------

- :doc:`API Reference <../../../api/core/split>` - Complete function signatures and detailed examples
- :doc:`Quartets <../quartets/overview>` - Quartet-based phylogenetic representations
- :doc:`Networks (Basic) <../networks/basic>` - Network classes using splits
