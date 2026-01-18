Overview
========

Quartets are four-taxon unrooted trees and are fundamental building blocks for 
network inference :cite:`PhyloZoo2024`. The quartet module provides classes for 
working with quartets, quartet profiles, and quartet profile sets.

Overview
--------

A quartet represents an unrooted tree on four taxa. There are three possible 
quartet topologies, or a star quartet (unresolved). Quartets are used extensively 
in network inference algorithms.

The quartet module provides three main classes:

* **Quartet**: Represents a single quartet (resolved or star)
* **QuartetProfile**: Represents a probability distribution over quartet topologies for a 4-taxon set
* **QuartetProfileSet**: A collection of quartet profiles with weights

Basic Usage
-----------

.. code-block:: python

   from phylozoo import Quartet, QuartetProfile, QuartetProfileSet
   from phylozoo.core.split import Split
   
   # Create resolved quartets
   q1 = Quartet(Split({"A", "B"}, {"C", "D"}))  # AB|CD
   q2 = Quartet(Split({"A", "C"}, {"B", "D"}))  # AC|BD
   
   # Create star quartet (unresolved)
   star = Quartet({"A", "B", "C", "D"})
   
   # Create profile
   profile = QuartetProfile(
       quartet1=q1,
       weight1=0.7,
       quartet2=q2,
       weight2=0.3
   )
   
   # Create profile set
   profile_set = QuartetProfileSet([profile])

For detailed information on each class, see:

* :doc:`Quartet <quartet>`: Single quartet operations
* :doc:`QuartetProfile <quartet_profile>`: Profile operations
* :doc:`QuartetProfileSet <quartet_profile_set>`: Profile set operations

.. seealso::
   For network inference using quartets, see :doc:`Inference <../../inference/overview>`. 
   For extracting quartets from networks, see :doc:`Networks (Advanced) <../networks/directed/advanced>`.
