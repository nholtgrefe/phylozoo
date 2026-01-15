Distance Matrices (Advanced)
=============================

This page covers advanced distance matrix operations including Kalmanson matrices, 
pseudo-metrics, and TSP solving. For basic operations, see 
:doc:`Distance (Basic) <core_distance_basic>`.

Advanced Classifications
-------------------------

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
   Kalmanson matrices are important for quartet-based network inference methods. If 
   your distance matrix is Kalmanson, it indicates a specific evolutionary structure.

.. tip::
   Use approximate TSP tour for large distance matrices. The nearest neighbor heuristic 
   is fast and often produces good results.

.. warning::
   Optimal TSP tour has exponential time complexity. Use only for small matrices 
   (< 15 taxa). For larger matrices, use approximate methods.
