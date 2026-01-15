Distance Matrices (Basic)
=========================

This page covers basic operations for working with distance matrices. For advanced 
classifications and TSP operations, see :doc:`Distance (Advanced) <core_distance_advanced>`.

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

See the :doc:`I/O <io>` page for supported formats.

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

.. note::
   Distance matrices are immutable. To modify a distance matrix, create a new instance 
   with the modified data.

.. tip::
   Use basic classifications to validate distance matrices before using them in 
   phylogenetic analyses. Many methods require metric or pseudo-metric distance matrices.

.. seealso::
   For advanced classifications (Kalmanson, pseudo-metric) and TSP operations, see 
   :doc:`Distance (Advanced) <core_distance_advanced>`.
