Primitives
==========

The primitives module (`phylozoo.core.primitives`) provides fundamental data structures 
used throughout PhyloZoo :cite:`PhyloZoo2024`. These classes form the building blocks 
for more complex phylogenetic data types.

Partition
---------

A ``Partition`` represents a partition of a set into disjoint subsets:

.. code-block:: python

   from phylozoo.core.primitives import Partition
   
   # Create a partition
   partition = Partition([{1, 2}, {3, 4}, {5}])
   
   # Access parts
   parts = partition.parts  # Frozen set of frozensets
   elements = partition.elements  # All elements in the partition
   
   # Check if partition
   is_valid = partition.is_partition()  # True
   
   # Check refinement
   other = Partition([{1}, {2}, {3, 4}, {5}])
   is_refinement = partition.is_refinement(other)  # True

Circular Orderings
------------------

Circular orderings represent elements arranged in a circle:

.. code-block:: python

   from phylozoo.core.primitives import CircularOrdering, CircularSetOrdering
   
   # Create circular ordering of elements
   co = CircularOrdering([1, 2, 3, 4])
   
   # Check if elements are neighbors
   are_neighbors = co.are_neighbors(1, 2)  # True
   are_neighbors = co.are_neighbors(1, 4)  # True (circular)
   
   # Create circular ordering of sets
   cso = CircularSetOrdering([{1, 2}, {3}, {4, 5}])
   
   # Check if sets are neighbors
   are_neighbors = cso.are_neighbors({1, 2}, {3})  # True

Circular orderings are stored in canonical form, making equality comparisons 
efficient. Two orderings are equal if they are the same up to cyclic permutation 
or reversal.

Graph Primitives
----------------

PhyloZoo provides graph primitives for working with network structures:

.. code-block:: python

   from phylozoo.core.primitives import DirectedMultiGraph, MixedMultiGraph
   
   # Create directed multigraph (supports parallel edges)
   dmg = DirectedMultiGraph(edges=[(1, 2), (1, 2)])  # Parallel edges
   
   # Create mixed multigraph (directed and undirected edges)
   mmg = MixedMultiGraph(
       edges=[(1, 2, "directed"), (2, 3, "undirected")]
   )

These graph primitives are used internally by network classes but can also be 
used directly for graph operations.

Available Classes
-----------------

**Partition:**

* **Partition** - Immutable partition class representing a partition of a set into 
  disjoint subsets. Supports refinement checking, size computation, and element 
  access. See :class:`phylozoo.core.primitives.Partition` for full API.

**Circular Orderings:**

* **CircularOrdering** - Immutable circular ordering of elements. Each element is 
  in its own singleton set. Stored in canonical form for efficient equality. See 
  :class:`phylozoo.core.primitives.CircularOrdering` for full API.

* **CircularSetOrdering** - Immutable circular ordering of sets. More general than 
  CircularOrdering, allows multiple elements per set. See 
  :class:`phylozoo.core.primitives.CircularSetOrdering` for full API.

**Graph Primitives:**

* **DirectedMultiGraph** - Directed multigraph supporting parallel edges. Used as 
  the base for DirectedPhyNetwork. See :class:`phylozoo.core.primitives.dmultigraph.DirectedMultiGraph` 
  for full API.

* **MixedMultiGraph** - Mixed multigraph supporting both directed and undirected 
  edges. Used as the base for SemiDirectedPhyNetwork. See 
  :class:`phylozoo.core.primitives.mmultigraph.MixedMultiGraph` for full API.

**Key Methods:**

* **Partition.is_refinement(other)** - Check if this partition is a refinement of 
  another partition. Returns boolean.

* **CircularOrdering.are_neighbors(a, b)** - Check if two elements are neighbors in 
  the circular ordering. Returns boolean.

* **CircularSetOrdering.are_neighbors(set1, set2)** - Check if two sets are neighbors 
  in the circular ordering. Returns boolean.

.. note::
   Primitives are typically used internally by higher-level classes. You may not 
   need to use them directly unless you're implementing custom functionality.

.. tip::
   Circular orderings are particularly useful for quartet-based methods and distance 
   matrix classifications (e.g., Kalmanson matrices).

.. seealso::
   For network classes that use these primitives, see :doc:`Networks (Basic) <core_networks_basic>`. 
   For split systems (which use Partition), see :doc:`Splits <core_splits>`.
