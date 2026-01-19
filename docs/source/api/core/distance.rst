Distance
========

The distance module provides classes and functions for working with distance matrices.
A distance matrix represents pairwise distances between a set of labeled items, where
distances satisfy properties such as symmetry :math:`d(i,j) = d(j,i)` and non-negativity
:math:`d(i,j) \geq 0`.

Main Classes
------------

.. automodule:: phylozoo.core.distance.base
   :members:
   :show-inheritance:

Classification Functions
------------------------

The following functions classify distance matrices based on mathematical properties:

- **Triangle inequality**: A distance matrix satisfies the triangle inequality if
  :math:`d(i,k) \leq d(i,j) + d(j,k)` for all :math:`i, j, k`.
- **Metric properties**: A metric is a distance function that satisfies the triangle
  inequality, symmetry, and non-negativity.
- **Kalmanson conditions**: A distance matrix is Kalmanson if it satisfies certain
  circular ordering constraints.

.. automodule:: phylozoo.core.distance.classifications
   :members:
   :show-inheritance:

Operations
----------

The operations module provides algorithms for working with distance matrices, including
the Traveling Salesman Problem (TSP) solver using the Held-Karp dynamic programming
algorithm.

.. automodule:: phylozoo.core.distance.operations
   :members:
   :show-inheritance:

I/O Support
-----------

Distance matrices support reading and writing in multiple formats: NEXUS, PHYLIP, and CSV.

.. automodule:: phylozoo.core.distance.io
   :members:
   :show-inheritance:
