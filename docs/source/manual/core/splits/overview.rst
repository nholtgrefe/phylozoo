Overview
========

Splits represent bipartitions of taxa and are a fundamental way to encode phylogenetic 
relationships :cite:`PhyloZoo2024`. The split module provides classes for working with 
splits and split systems.

Overview
--------

A split represents a bipartition of a set of elements. Splits are used to encode 
phylogenetic relationships and can be extracted from networks or used to construct trees.

The split module provides two main classes:

* **Split**: Represents a single bipartition
* **SplitSystem**: A collection of splits covering a complete set of elements

Basic Usage
-----------

.. code-block:: python

   from phylozoo import Split, SplitSystem
   
   # Create a split
   split = Split({"A", "B"}, {"C", "D"})
   
   # Check compatibility
   other_split = Split({"A"}, {"B", "C", "D"})
   are_compatible = split.is_compatible(other_split)
   
   # Create split system
   splits = [
       Split({"A", "B"}, {"C", "D"}),
       Split({"A"}, {"B", "C", "D"}),
   ]
   split_system = SplitSystem(splits)

For detailed information on each class, see:

* :doc:`Split <split>`: Single split operations
* :doc:`SplitSystem <split_system>`: Split system operations

.. seealso::
   For extracting splits from networks, see :doc:`Networks (Advanced) <../networks/directed/advanced>`. 
   For I/O operations, see :doc:`I/O <../../io>`.
