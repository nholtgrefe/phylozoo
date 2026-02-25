PHYLIP Format
=============

PHYLIP is a compact format for distance matrices, compatible with standard PHYLIP
tools. It is space-efficient and widely used for phylogenetic distance data.

**Classes using PHYLIP format:** :class:`phylozoo.core.distance.DistanceMatrix`

**File extensions:** ``.phy``, ``.phylip``

Format Structure
----------------

The first line is the number of taxa. Each following line contains a taxon name
(padded to 10 characters in strict PHYLIP) followed by the distances for that taxon:

.. code-block:: text

   3
   A          0.00000 1.00000 2.00000
   B          1.00000 0.00000 1.00000
   C          2.00000 1.00000 0.00000

Example
-------

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np

   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   dm.save("distances.phy", format="phylip")
   dm2 = DistanceMatrix.load("distances.phy", format="phylip")

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`nexus` — NEXUS distance matrix format
   * :doc:`csv` — Tabular distance matrix format
   * :doc:`../../../core/distance` — Distance matrices
