CSV Format
=========

CSV (comma-separated values) is a simple tabular format for distance matrices, easy
to read and edit in spreadsheet applications. PhyloZoo supports a header row with
taxon labels and optional custom delimiters (e.g. tab for TSV).

**Classes using CSV format:** :class:`phylozoo.core.distance.DistanceMatrix`

**File extensions:** ``.csv`` (and optionally ``.tsv`` with ``delimiter="\t"``)

Format Structure
----------------

The first row is a header: an empty first cell followed by taxon labels. Each
following row has a taxon label in the first column and then the distances:

.. code-block:: text

   ,A,B,C
   A,0.0,1.0,2.0
   B,1.0,0.0,1.0
   C,2.0,1.0,0.0

Example
-------

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np

   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   dm.save("distances.csv", format="csv")
   dm2 = DistanceMatrix.load("distances.csv", format="csv")
   dm.save("distances.tsv", format="csv", delimiter="\t")

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`nexus` — NEXUS distance format
   * :doc:`phylip` — PHYLIP distance format
   * :doc:`../../../core/distance` — Distance matrices
