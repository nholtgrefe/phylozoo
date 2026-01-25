Overview
========

The :mod:`phylozoo.core.primitives` module provides fundamental, immutable data structures
that serve as building blocks for phylogenetic analysis in PhyloZoo. These classes
represent core mathematical concepts used throughout the library, including partitions,
circular orderings, and graph primitives that form the foundation for network-based
phylogenetics.

All classes and functions on this page can be imported from the core primitives module:

.. code-block:: python

   from phylozoo.core.primitives import *
   # or directly
   from phylozoo.core.primitives import Partition, CircularOrdering

Overview of Primitive Classes
-----------------------------

This section provides brief descriptions of each primitive data structure. For detailed
documentation and usage examples, see the individual pages linked below.

**Partition** - Represents a partition of a set into disjoint subsets. Used as the base
for split systems and provides methods for checking refinement relationships and other
partition operations. Essential for analyzing phylogenetic splits and compatibility.

**Circular Ordering** - Represents elements arranged in a circle, with each element in
its own singleton set. Stored in canonical form for efficient equality comparisons.
Fundamental for quartet-based methods and distance matrix classifications.

**Circular Set Ordering** - A generalization of circular ordering where elements can be
grouped into sets. More flexible than CircularOrdering and used for complex phylogenetic
arrangements where elements naturally cluster together.

**Directed Multi-Graph** - A directed multigraph that supports parallel edges between
nodes. Serves as the foundation for DirectedPhyNetwork and enables representation of
complex directed phylogenetic relationships with reticulations.

**Mixed Multi-Graph** - A mixed multigraph supporting both directed and undirected edges,
plus parallel edges of both types. Forms the basis for SemiDirectedPhyNetwork and allows
modeling phylogenetic networks with both tree-like and reticulate relationships.

.. note::
   :class: dropdown

   **Implementation details**

   Primitive data structures are designed for immutability and mathematical precision:

   - All primitives are immutable after construction (frozen sets, tuples, canonical forms)
   - Circular orderings use canonical representations for efficient equality testing
   - Partitions maintain mathematical properties (disjoint, covering) through validation
   - Graph primitives extend NetworkX with phylogenetic-specific operations
   - Memory-efficient representations minimize overhead for large phylogenetic datasets

   For implementation details, see :mod:`src/phylozoo/core/primitives/`.

Working with Primitives
------------------------

Primitives are typically used internally by higher-level PhyloZoo classes, but can also
be used directly for custom phylogenetic analyses. They provide mathematically rigorous
foundations for:

- **Split decomposition** and **phylogenetic networks** (Partition)
- **Quartet systems** and **distance classifications** (Circular Ordering)
- **Network construction** and **reticulate evolution** (Graph primitives)

Most users will interact with primitives indirectly through network classes, but direct
use is recommended when implementing custom phylogenetic algorithms or analyzing
mathematical properties of phylogenetic data.

.. tip::
   Circular orderings are particularly useful for quartet-based methods and distance
   matrix classifications (e.g., Kalmanson matrices). Partitions excel at analyzing
   phylogenetic splits and compatibility.

See Also
--------

- :doc:`Partition <partition>` - Working with set partitions and refinement relationships
- :doc:`Circular Ordering <circular_ordering>` - Circular arrangements of elements and sets
- :doc:`Directed Multi-Graph <directed_multigraph>` - Directed graphs with parallel edges
- :doc:`Mixed Multi-Graph <mixed_multigraph>` - Mixed graphs for phylogenetic networks
- :doc:`API Reference <../../../api/core/primitives>` - Complete function signatures and detailed examples
- :doc:`Networks (Basic) <../networks/basic>` - Network classes that use these primitives
- :doc:`Splits <../splits/overview>` - Split systems based on partitions
