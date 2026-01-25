Pairwise Distances
==================

The :mod:`phylozoo.core.distance` module provides immutable, well-typed containers for pairwise distances
between taxa, along with comprehensive tools for classification and optimization.
Distance matrices are fundamental to many phylogenetic inference algorithms, serving
as the standard input format for quartet profile generation, network reconstruction
algorithms such as Squirrel, and other downstream analyses.


All classes and functions on this page can be imported from the core distance module:

.. code-block:: python

   from phylozoo.core.distance import *
   # or directly
   from phylozoo.core.distance import DistanceMatrix

Working with Distance Matrices
--------------------------------

The :class:`phylozoo.core.distance.DistanceMatrix` class is the canonical container
for pairwise distances in PhyloZoo. It provides an immutable, labeled, and read-only
representation that ensures data integrity throughout your analysis pipeline.

.. note::
   :class: dropdown

   **Implementation details**

   Distance matrices are designed for immutability and performance:

   - The underlying NumPy array is made read-only after construction
   - Labels are stored as an immutable tuple and must be hashable (e.g., ``str``, ``int``, ``tuple``)
   - To modify distances, create a new ``DistanceMatrix`` instance with the updated data
   - The matrix is validated for symmetry and squareness at construction time

   For implementation details, see :mod:`src/phylozoo/core/distance/base.py`.

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
and raises :class:`phylozoo.utils.exceptions.PhyloZooValueError` if the input is invalid.

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

- **NEXUS** (default): Standard phylogenetic data format
- **PHYLIP**: Classic phylogenetic format
- **CSV**: Comma-separated values for interoperability

.. code-block:: python

   # Load from file (auto-detects format by extension)
   dm = DistanceMatrix.load("distances.nexus")

   # Load with explicit format
   dm = DistanceMatrix.load("distances.phy", format="phylip")

   # Save to file
   dm.save("output.csv", format="csv")

.. seealso::
   The I/O system uses the :class:`phylozoo.utils.io.IOMixin` interface, providing
   consistent file handling across PhyloZoo classes. For details on the I/O system,
   see the :doc:`I/O documentation <../utils/io>`. For specific information about
   supported file formats and parameter options for distance matrices, see the
   :mod:`API reference <phylozoo.core.distance.io>`.

Classification Functions
-------------------------

The distance module provides functions to classify distance matrices based on
mathematical properties. These classifications are essential for validating
that distance matrices meet the requirements of specific algorithms.

Basic Properties
^^^^^^^^^^^^^^^^

**Non-negativity**

The :func:`phylozoo.core.distance.classifications.is_nonnegative` function checks
that all distances are non-negative: :math:`d(i,j) \geq 0` for all :math:`i, j`.
This is a fundamental requirement for distance functions.

**Zero Diagonal**

The :func:`phylozoo.core.distance.classifications.has_zero_diagonal` function
verifies that all diagonal entries are exactly zero, ensuring that the distance
from any item to itself is zero. This is a strict equality check using NumPy
equality, intended to verify metric conventions.

**Triangle Inequality**

The :func:`phylozoo.core.distance.classifications.satisfies_triangle_inequality`
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

The :func:`phylozoo.core.distance.classifications.is_metric` function provides
a convenient wrapper that checks all of these properties. Use this function when
downstream algorithms require a proper metric and you want a single boolean check.

**Pseudo-metrics**

A pseudo-metric is like a metric but allows distinct items to have zero distance.
The :func:`phylozoo.core.distance.classifications.is_pseudo_metric` function
checks non-negativity and triangle inequality, but does not require a zero diagonal.
This is useful when dealing with distance matrices that may have identical items.

**Kalmanson Conditions**

The Kalmanson condition is a pair of inequalities that must hold for every ordered
quadruple :math:`(i < j < k < l)` in a circular ordering:

.. math::

   d(i,j) + d(k,l) &\leq d(i,k) + d(j,l) \\
   d(i,l) + d(j,k) &\leq d(i,k) + d(j,l)

The :func:`phylozoo.core.distance.classifications.is_kalmanson` function checks
these conditions with respect to a given circular ordering. This function requires:

- The matrix to be a pseudo-metric
- The provided :class:`phylozoo.core.primitives.CircularOrdering` to contain
  exactly the same labels as the distance matrix

The check is implemented using Numba-accelerated loops and will raise
:class:`phylozoo.utils.exceptions.PhyloZooValueError` for invalid inputs.

Traveling Salesman Problem (TSP)
---------------------------------

The distance module provides functions for solving the Traveling Salesman Problem,
which finds the optimal tour visiting all taxa exactly once and returning to the
starting point. TSP solutions are useful for generating circular orderings and
optimizing distance-based analyses.

Exact Solution
^^^^^^^^^^^^^^

The :func:`phylozoo.core.distance.operations.optimal_tsp_tour` function provides
an exact solution using the Held–Karp dynamic programming algorithm :cite:`HeldKarp1962`.

The function returns a :class:`phylozoo.core.primitives.CircularOrdering` in canonical
form, representing the optimal tour.

Example:

.. code-block:: python

   from phylozoo.core.distance.operations import optimal_tsp_tour

   tour = optimal_tsp_tour(dm)
   print(f"Optimal tour: {tour.order}")

Use this function when exact optimality is required and the problem size is manageable.

.. note::
   :class: dropdown

   **Implementation details**

   The Held–Karp algorithm implementation:

   - Uses bitmask-based dynamic programming for efficient state representation
   - Is accelerated with Numba JIT compilation for performance
   - Memory grows exponentially with :math:`2^n`, so avoid large :math:`n`
   - For large instances, consider using approximate methods instead

   The algorithm uses a memoization table of size :math:`n \times 2^{n-1}` to store
   optimal subproblem solutions, making it practical only for small instances
   (typically :math:`n \leq 15`–:math:`20` depending on available memory).

Approximate Solutions
^^^^^^^^^^^^^^^^^^^^^

For larger instances where exact solutions are infeasible, the
:func:`phylozoo.core.distance.operations.approximate_tsp_tour` function provides
heuristic solvers. Supported methods include:

- **``simulated_annealing``** (default): Simulated annealing with greedy initialization.
  Often produces good-quality tours at moderate runtime.

- **``greedy``**: Nearest-neighbor heuristic. Very fast but can produce poor results
  in adversarial cases.

- **``christofides``**: Christofides algorithm providing a :math:`3/2`-approximation
  when the distance matrix is metric. Implemented using NetworkX utilities.

The function uses NetworkX's traveling salesman utilities and returns a
:class:`phylozoo.core.primitives.CircularOrdering`. It is suitable for larger
matrices where exact solutions are computationally prohibitive.

See Also
--------

- :doc:`API Reference <../../api/core/distance>` - Complete function signatures and detailed examples
- :mod:`phylozoo.core.distance.io` - Distance matrix I/O format details and parameter options
