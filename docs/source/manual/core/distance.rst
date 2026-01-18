Distance Matrices
=================

This page covers working with distance matrices in PhyloZoo, including basic operations, 
classifications, and TSP operations.

**I/O Formats**: NEXUS (default, extensions: ``.nexus``, ``.nex``, ``.nxs``), PHYLIP (extensions: ``.phy``, ``.phylip``), 
CSV (extension: ``.csv``). See :doc:`I/O <../io>` for details.

Creating Distance Matrices
--------------------------

Create distance matrices from NumPy arrays:

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np
   
   # Create from numpy array
   matrix = np.array([
       [0, 1, 2, 3],
       [1, 0, 1, 2],
       [2, 1, 0, 1],
       [3, 2, 1, 0]
   ])
   
   dm = DistanceMatrix(matrix, labels=["A", "B", "C", "D"])

Distance matrices can also be loaded from files:

.. code-block:: python

   # Load from file
   dm = DistanceMatrix.load("distances.nexus")
   dm = DistanceMatrix.load("distances.phy", format="phylip")

See the :doc:`I/O <../reference/io>` page for supported formats.

Accessing Distances
-------------------

Access distances between taxa:

.. code-block:: python

   # Get distance between two taxa
   dist = dm.get_distance("A", "B")
   
   # Get index for a label
   idx = dm.get_index("A")
   
   # Access underlying numpy array
   array = dm.np_array

Basic Classifications
---------------------

Check basic distance matrix properties:

.. code-block:: python

   from phylozoo.core.distance import classifications
   
   # Check triangle inequality
   is_metric = classifications.satisfies_triangle_inequality(dm)
   
   # Check zero diagonal
   has_zero = classifications.has_zero_diagonal(dm)
   
   # Check symmetry
   is_symmetric = classifications.is_symmetric(dm)

Advanced Classifications
------------------------

Classify distance matrices for specific phylogenetic properties:

.. code-block:: python

   from phylozoo.core.distance import classifications
   from phylozoo.core.primitives import CircularOrdering
   
   # Check if pseudo-metric
   is_pseudometric = classifications.is_pseudo_metric(dm)
   
   # Check if Kalmanson (circular decomposable)
   ordering = CircularOrdering(["A", "B", "C", "D"])
   is_kalmanson = classifications.is_kalmanson(dm, ordering)
   
   # Check if ultrametric
   is_ultrametric = classifications.is_ultrametric(dm)

Kalmanson matrices are particularly important for certain phylogenetic methods as they 
indicate that the distance matrix can be decomposed along a circular ordering.

TSP Operations
--------------

Solve Traveling Salesman Problem (TSP) to find optimal circular orderings:

.. code-block:: python

   from phylozoo.core.distance import operations
   from phylozoo.core.primitives import CircularOrdering
   
   # Find optimal TSP tour (exact solution, may be slow for large matrices)
   optimal_ordering = operations.optimal_tsp_tour(dm)
   
   # Approximate TSP tour (faster, good for large matrices)
   approx_ordering = operations.approximate_tsp_tour(dm, method="nearest_neighbor")

TSP tours are useful for finding circular orderings that minimize total distance, 
which can be used with Kalmanson classification methods.

Available Functions
-------------------

**Classes:**

* **DistanceMatrix** - Immutable distance matrix class. Represents pairwise distances 
  between labeled items. Stored as symmetric NumPy array. Supports I/O operations and 
  various classifications. See :class:`phylozoo.core.distance.DistanceMatrix` for 
  full API.

**Basic Methods:**

* **get_distance(label1, label2)** - Get distance between two labels. Returns float.
* **get_index(label)** - Get index for a label. Returns integer.
* **labels** - Tuple of labels (immutable).
* **np_array** - Read-only access to underlying NumPy array.
* **save(filename)** - Save distance matrix to file.
* **load(filename)** - Load distance matrix from file (class method).

**Basic Classifications** (``phylozoo.core.distance.classifications``):

* **satisfies_triangle_inequality(distance_matrix)** - Check if distance matrix 
  satisfies triangle inequality (d(i,k) <= d(i,j) + d(j,k) for all i, j, k). Returns 
  boolean. See :func:`phylozoo.core.distance.classifications.satisfies_triangle_inequality`.

* **has_zero_diagonal(distance_matrix)** - Check if diagonal is zero. Returns boolean.

* **is_symmetric(distance_matrix)** - Check if matrix is symmetric. Returns boolean.

**Advanced Classifications** (``phylozoo.core.distance.classifications``):

* **is_pseudo_metric(distance_matrix)** - Check if distance matrix is a pseudo-metric. 
  A pseudo-metric satisfies triangle inequality and has zero diagonal, but distances 
  can be zero for distinct items. Returns boolean. See 
  :func:`phylozoo.core.distance.classifications.is_pseudo_metric`.

* **is_kalmanson(distance_matrix, ordering)** - Check if distance matrix is Kalmanson 
  with respect to a circular ordering. Kalmanson matrices satisfy specific inequalities 
  that make them decomposable along the circular order. Returns boolean. See 
  :func:`phylozoo.core.distance.classifications.is_kalmanson`.

* **is_ultrametric(distance_matrix)** - Check if distance matrix is ultrametric. An 
  ultrametric satisfies the strong triangle inequality: d(i,k) <= max(d(i,j), d(j,k)) 
  for all i, j, k. Returns boolean. See 
  :func:`phylozoo.core.distance.classifications.is_ultrametric`.

**TSP Operations** (``phylozoo.core.distance.operations``):

* **optimal_tsp_tour(distance_matrix)** - Find optimal TSP tour using Held-Karp 
  algorithm. Returns CircularOrdering. Exact solution but exponential time complexity. 
  Use for small matrices (typically < 15 taxa). See 
  :func:`phylozoo.core.distance.operations.optimal_tsp_tour`.

* **approximate_tsp_tour(distance_matrix, method="nearest_neighbor")** - Find 
  approximate TSP tour using heuristic methods. Returns CircularOrdering. Much faster 
  than optimal solution. Methods: "nearest_neighbor", "2opt". See 
  :func:`phylozoo.core.distance.operations.approximate_tsp_tour`.

.. note::
   Distance matrices are immutable. To modify a distance matrix, create a new instance 
   with the modified data.

.. note::
   Kalmanson matrices are important for quartet-based network inference methods. If 
   your distance matrix is Kalmanson, it indicates a specific evolutionary structure.

.. tip::
   Use basic classifications to validate distance matrices before using them in 
   phylogenetic analyses. Many methods require metric or pseudo-metric distance matrices.

.. tip::
   Use approximate TSP tour for large distance matrices. The nearest neighbor heuristic 
   is fast and often produces good results.

.. warning::
   Optimal TSP tour has exponential time complexity. Use only for small matrices 
   (< 15 taxa). For larger matrices, use approximate methods.
