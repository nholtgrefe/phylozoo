CSV
===

CSV (comma-separated values) is a simple tabular format for distance matrices, easy
to read and edit in spreadsheet applications. PhyloZoo supports a header row with
taxon labels and optional custom delimiters (e.g. tab for TSV).

.. seealso::
   `Comma-separated values <https://en.wikipedia.org/wiki/Comma-separated_values>`_ — Wikipedia

Classes and extensions
----------------------

**Classes:** :class:`~phylozoo.core.distance.base.DistanceMatrix`

**File extensions:** ``.csv`` (and optionally ``.tsv`` with ``delimiter="\t"``)

Structure
---------

The first row is a header: an empty first cell followed by taxon labels. Each
following row has a taxon label in the first column and then the distances:

.. code-block:: text

   ,A,B,C
   A,0.0,1.0,2.0
   B,1.0,0.0,1.0
   C,2.0,1.0,0.0

Examples
--------

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np

   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   dm.save("distances.csv", format="csv")
   dm2 = DistanceMatrix.load("distances.csv", format="csv")
   dm.save("distances.tsv", format="csv", delimiter="\t")

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`nexus` — NEXUS distance format
- :doc:`phylip` — PHYLIP distance format
- :doc:`../../../core/distance` — Distance matrices
