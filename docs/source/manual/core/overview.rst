Overview
========

The core module (`phylozoo.core`) contains fundamental data structures and classes used 
throughout PhyloZoo :cite:`PhyloZoo2024`. All core classes are immutable and designed for 
type safety and performance.

Overview
--------

The core module provides:

* **Networks**: Directed and semi-directed phylogenetic network classes
* **Quartets**: Four-taxon relationships and profiles
* **Splits**: Bipartitions and split systems
* **Sequences**: Multiple sequence alignment (MSA) handling
* **Distance**: Distance matrix classes and operations
* **Primitives**: Fundamental structures (Partition, CircularOrdering, etc.)

All core data structures are available at the package top-level for convenience:

.. code-block:: python

   from phylozoo import (
       DirectedPhyNetwork,
       SemiDirectedPhyNetwork,
       Quartet,
       QuartetProfile,
       QuartetProfileSet,
       Split,
       SplitSystem,
       DistanceMatrix,
       MSA,
   )

Key Components
--------------

The `phylozoo.core` module is organized into several submodules, each focusing on a 
specific data type or functionality:

**Networks:**
*   :doc:`Networks Overview <networks/overview>`: Overview of network classes
*   :doc:`Directed Networks <networks/directed/overview>`: Creating directed networks, basic properties,
    loading/saving, advanced features, generators
*   :doc:`Semi-Directed Networks <networks/semi_directed/overview>`: Creating semi-directed networks, basic properties,
    loading/saving, advanced features, generators

**Sequences:**
*   :doc:`Sequences <sequences>`: Creating MSAs, basic operations, bootstrapping,
    per-gene bootstrapping, efficient array operations

**Distance Matrices:**
*   :doc:`Distance Matrices <distance>`: Creating distance matrices, basic and advanced
    classifications, TSP operations

**Other Components:**
*   :doc:`Quartets Overview <quartets/overview>`: Four-taxon unrooted trees and profiles
*   :doc:`Splits Overview <splits/overview>`: Bipartitions of taxa and split systems
*   :doc:`Primitives Overview <primitives/overview>`: Fundamental data structures (Partition, CircularOrdering, etc.)

These components are designed to be interoperable, allowing for seamless transitions 
between different data representations.

I/O Support
-----------

All core classes support reading and writing to various file formats. See the 
:doc:`I/O <../io>` page for details on supported formats and usage. Common formats include:

* **eNewick**: Extended Newick format for networks (see :cite:`Cardona2008`)
* **FASTA/NEXUS**: For sequence alignments
* **NEXUS/PHYLIP/CSV**: For distance matrices
* **DOT**: Graphviz format for networks

Design Principles
-----------------

**Immutability**: All core classes are immutable after initialization. To modify a 
network or other structure, create a new instance with the desired changes.

**Type Safety**: Comprehensive type hints enable static type checking and better IDE 
support.

**Validation**: Objects are validated upon creation to ensure they represent valid 
phylogenetic structures.

For more detailed information on specific components, see the subpages above.
