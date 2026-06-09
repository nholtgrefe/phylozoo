Pairwise Distances
==================

The :mod:`phylozoo.core.distance` module provides immutable containers for pairwise distances
between taxa, along with comprehensive tools for classification and optimization.
Distance matrices are fundamental to many phylogenetic inference algorithms.


All classes and functions on this page can be imported from the core distance module:

.. code-block:: python

   from phylozoo.core.distance import *
   # or directly
   from phylozoo.core.distance import DistanceMatrix

Working with Distance Matrices
--------------------------------

The :class:`~phylozoo.core.distance.base.DistanceMatrix` class is the canonical container
for pairwise distances in PhyloZoo. It provides an immutable, labeled, and read-only
representation that ensures data integrity throughout your analysis pipeline.

Creating a Distance Matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Distance matrices can be created from NumPy arrays with optional labels:

.. code-block:: python

   import numpy as np
   from phylozoo.core.distance import DistanceMatrix

   # Create a symmetric distance matrix
   matrix = np.array([
       [0.0, 1.0, 2.0],
       [1.0, 0.0, 1.0],
       [2.0, 1.0, 0.0]
   ])

   # Create with labels
   dm = DistanceMatrix(matrix, labels=["A", "B", "C"])

   # Or without labels (defaults to 0, 1, 2, ...)
   dm2 = DistanceMatrix(matrix)

The matrix must be square and symmetric. The constructor validates these properties
and raises :class:`~phylozoo.utils.exceptions.general.PhyloZooValueError` if the input is invalid.

Accessing Distances
^^^^^^^^^^^^^^^^^^^

Distance matrices provide several methods for accessing distances and metadata:

.. code-block:: python

   # Get distance between two labeled items
   distance = dm.get_distance("A", "B")  # Returns: 1.0

   # Get the index of a label
   index = dm.get_index("A")  # Returns: 0

   # Access labels and underlying array
   labels = dm.labels        # Returns: ("A", "B", "C")
   array = dm.np_array       # Returns: read-only numpy array

These methods provide safe access to distances without exposing the mutable underlying
array, maintaining immutability guarantees.

File Input/Output
^^^^^^^^^^^^^^^^^

Distance matrices support reading and writing in multiple phylogenetic formats:

- **NEXUS** (default): Standard phylogenetic data format — see :doc:`NEXUS format <../utils/io/formats/nexus>`
- **PHYLIP**: Classic phylogenetic format — see :doc:`PHYLIP format <../utils/io/formats/phylip>`
- **CSV**: Comma-separated values for interoperability — see :doc:`CSV format <../utils/io/formats/csv>`

.. code-block:: python

   # Load from file (auto-detects format by extension)
   dm = DistanceMatrix.load("distances.nexus")

   # Load with explicit format
   dm = DistanceMatrix.load("distances.phy", format="phylip")

   # Save to file
   dm.save("output.csv", format="csv")

.. seealso::   
   The `DistanceMatrix` class uses the :class:`~phylozoo.utils.io.mixin.IOMixin` interface, providing
   consistent file handling across PhyloZoo classes. For details on the I/O system,
   see the :doc:`I/O manual <../utils/io/overview>`.

Classification Functions
-------------------------

The distance module provides functions to classify distance matrices based on
mathematical properties. These classifications are essential for validating
that distance matrices meet the requirements of specific algorithms.

Basic Properties
^^^^^^^^^^^^^^^^

**Non-negativity**

The :func:`~phylozoo.core.distance.classifications.is_nonnegative` function checks
that all distances are non-negative: :math:`d(i,j) \geq 0` for all :math:`i, j`.
This is a fundamental requirement for distance functions.

**Zero Diagonal**

The :func:`~phylozoo.core.distance.classifications.has_zero_diagonal` function
verifies that all diagonal entries are exactly zero, ensuring that the distance
from any item to itself is zero. This is a strict equality check using NumPy
equality, intended to verify metric conventions.

**Triangle Inequality**

The :func:`~phylozoo.core.distance.classifications.satisfies_triangle_inequality`
function checks whether the distance matrix satisfies the triangle inequality:

.. math::

   d(i,k) \leq d(i,j) + d(j,k) \quad \forall i, j, k

This property is fundamental to metric spaces and is required by many distance-based
algorithms. The implementation uses Numba-accelerated triple loops for efficient
computation on moderately sized matrices.

Example:

.. code-block:: python

   from phylozoo.core.distance.classifications import satisfies_triangle_inequality

   if satisfies_triangle_inequality(dm):
       print("Matrix satisfies triangle inequality")
   else:
       print("Matrix violates triangle inequality")

Metric Properties
^^^^^^^^^^^^^^^^^

**Metrics**

A metric distance matrix satisfies four properties:

1. Non-negativity: :math:`d(x, y) \geq 0` for all :math:`x, y`
2. Triangle inequality: :math:`d(x, z) \leq d(x, y) + d(y, z)` for all :math:`x, y, z`
3. Zero diagonal: :math:`d(x, x) = 0` for all :math:`x`
4. Symmetry: :math:`d(x, y) = d(y, x)` for all :math:`x, y` (enforced by constructor)

The :func:`~phylozoo.core.distance.classifications.is_metric` function provides
a convenient wrapper that checks all of these properties. Use this function when
downstream algorithms require a proper metric and you want a single boolean check.

**Pseudo-metrics**

A pseudo-metric is like a metric but allows distinct items to have zero distance.
The :func:`~phylozoo.core.distance.classifications.is_pseudo_metric` function
checks non-negativity and triangle inequality, but does not require a zero diagonal.

**Kalmanson Conditions**

The Kalmanson condition :cite:`Kalmanson1975` is a pair of inequalities that must hold for every ordered
quadruple :math:`(i < j < k < l)` in a circular ordering:

.. math::

   d(i,j) + d(k,l) &\leq d(i,k) + d(j,l) \\
   d(i,l) + d(j,k) &\leq d(i,k) + d(j,l)

The :func:`~phylozoo.core.distance.classifications.is_kalmanson` function checks
these conditions with respect to a given circular ordering. This function requires:

- The matrix to be a pseudo-metric
- The provided :class:`~phylozoo.core.primitives.circular_ordering.CircularOrdering` to contain
  exactly the same labels as the distance matrix

The check is implemented using Numba-accelerated loops.

**Tree Metrics**

The :func:`~phylozoo.core.distance.classifications.is_tree_metric` function checks whether
a distance matrix is a *tree metric* — a metric that arises as the path-length distances on
some edge-weighted tree. The check is based on the *four-point condition*: for every choice
of four taxa the maximum of the three pairwise sums

.. math::

   \{d_{ij}+d_{kl},\; d_{ik}+d_{jl},\; d_{il}+d_{jk}\}

must be achieved by at least two of them :cite:`Bandelt1992`. The test uses a
Numba-accelerated loop over all :math:`\binom{n}{4}` quartets and returns in
:math:`O(n^4)` time.

.. code-block:: python

   from phylozoo.core.distance.classifications import is_tree_metric

   if is_tree_metric(dm):
       print("Distance matrix is a tree metric")

**Totally Decomposable Metrics**

The :func:`~phylozoo.core.distance.classifications.is_totally_decomposable` function checks
whether a distance matrix is *totally decomposable*: the split-prime residual :math:`d^0`
of the canonical split decomposition is zero, meaning the metric is expressed exactly as a
weighted sum of split metrics.  Every tree metric is totally decomposable, but the converse
does not hold — totally decomposable metrics may have incompatible d-splits. See next
section for details on split decompositions.

.. code-block:: python

   from phylozoo.core.distance.classifications import is_totally_decomposable

   if is_totally_decomposable(dm):
       print("All pairwise distances are explained by d-splits alone")


Split Decomposition
-------------------

The :mod:`phylozoo.core.distance.decomposition` module implements the *canonical split
decomposition* of :cite:`Bandelt1992`.  Every finite metric :math:`d` decomposes
uniquely as

.. math::

   d = d^0 + \sum_{S} \alpha_S\,\delta_S

where :math:`S` ranges over all *d-splits* (bipartitions with positive isolation index),
:math:`\alpha_S` is the isolation index of :math:`S`, :math:`\delta_S` is the
corresponding split metric, and :math:`d^0` is the *split-prime residual* — a metric that
admits no further d-splits.

All classes and functions can be imported from the core distance module:

.. code-block:: python

   from phylozoo.core.distance import isolation_index, split_decomposition

Isolation Index
^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.distance.decomposition.isolation_index` function computes the
isolation index of a split :math:`(A, B)`:

.. math::

   \alpha_{A,B} = \tfrac{1}{2}\,
       \min_{\substack{i,j\in A \\ k,l\in B}}
       \bigl(\max\{d_{ij}+d_{kl},\, d_{ik}+d_{jl},\, d_{il}+d_{jk}\} - d_{ij} - d_{kl}\bigr)

The index is always non-negative; it is strictly positive if and only if the split is
a *d-split* of the distance matrix.

.. code-block:: python

   import numpy as np
   from phylozoo.core.distance import DistanceMatrix, isolation_index
   from phylozoo.core.split import Split

   d = np.array([
       [ 0,  4,  5,  7, 13,  8,  6],
       [ 4,  0,  1,  3,  9, 12, 10],
       [ 5,  1,  0,  2,  8, 13, 11],
       [ 7,  3,  2,  0,  6, 11, 13],
       [13,  9,  8,  6,  0,  5,  7],
       [ 8, 12, 13, 11,  5,  0,  2],
       [ 6, 10, 11, 13,  7,  2,  0],
   ], dtype=float)
   dm = DistanceMatrix(d, labels=list("ABCDEFG"))

   s = Split({"E", "F", "G"}, {"A", "B", "C", "D"})
   print(isolation_index(dm, s))  # 6.0

Split Decomposition
^^^^^^^^^^^^^^^^^^^

The :func:`~phylozoo.core.distance.decomposition.split_decomposition` function computes the
full decomposition and returns both the :class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem`
of all d-splits and the residual :class:`~phylozoo.core.distance.base.DistanceMatrix`.

.. code-block:: python

   system, residual = split_decomposition(dm)

   print(len(system.splits))          # 4  (Table 1 has four d-splits)
   import numpy as np
   print(np.allclose(residual.np_array, 0))  # True  (totally decomposable)

   for sp in system.splits:
       print(sp, "→ α =", system.get_weight(sp))

Round-Trip Identity
^^^^^^^^^^^^^^^^^^^

The decomposition satisfies :math:`d^1 + d^0 = d` exactly (up to floating-point noise),
where :math:`d^1 = \sum_S \alpha_S\,\delta_S` is computed from the weighted split system
via :func:`~phylozoo.core.split.algorithms.distances_from_splitsystem`:

.. code-block:: python

   from phylozoo.core.split.algorithms import distances_from_splitsystem

   d1_dm = distances_from_splitsystem(system)
   labels = list(dm.labels)
   n = len(labels)
   d1 = np.array([[d1_dm.get_distance(labels[r], labels[c]) for c in range(n)]
                  for r in range(n)])

   assert np.allclose(d1 + residual.np_array, dm.np_array)  # always holds

Conversely, when a :class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem`
encodes a tree metric (all d-splits are pairwise compatible), the decomposition is the
inverse of :func:`~phylozoo.core.split.algorithms.distances_from_splitsystem`:

.. code-block:: python

   from phylozoo.core.split import WeightedSplitSystem, Split
   from phylozoo.core.distance import split_decomposition
   from phylozoo.core.split.algorithms import distances_from_splitsystem

   # Build a tree split system, convert to distances, then recover the system
   tree_system = WeightedSplitSystem({
       Split({"A"}, {"B", "C", "D"}): 1.0,
       Split({"B"}, {"A", "C", "D"}): 1.0,
       Split({"C"}, {"A", "B", "D"}): 1.0,
       Split({"D"}, {"A", "B", "C"}): 1.0,
       Split({"A", "B"}, {"C", "D"}): 2.0,
   })
   dm = distances_from_splitsystem(tree_system)
   recovered, residual = split_decomposition(dm)

   assert recovered.splits == tree_system.splits        # same d-splits
   assert np.allclose(residual.np_array, 0)             # zero residual

See Also
--------

- :doc:`API Reference <../../api/core/distance>` — Complete function signatures and detailed examples
- :doc:`Split Systems <splits/split_system>` — Working with weighted split systems and round-trips
- :doc:`Circular Orderings <primitives/circular_ordering>` — Working with circular orderings