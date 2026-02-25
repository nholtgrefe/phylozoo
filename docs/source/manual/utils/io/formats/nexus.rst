NEXUS Format
============

The NEXUS format is a flexible, block-structured format used widely in phylogenetics.
A NEXUS file starts with ``#NEXUS`` and contains one or more blocks, each of the form
``BEGIN blockname; ... END;``. PhyloZoo uses different block types for different
data kinds (distance matrices, sequence alignments, split systems). The same file
extensions are used for all NEXUS subtypes; the class you load or save with
determines which block structure is expected.

**File extensions:** ``.nexus``, ``.nex``, ``.nxs``

General Structure
-----------------

All NEXUS files begin with the ``#NEXUS`` token. Data are organized in blocks:

.. code-block:: text

   #NEXUS

   BEGIN TAXA;
       ...
   END;

   BEGIN SomeBlock;
       ...
   END;

PhyloZoo supports the following NEXUS subtypes (block types).

Distance Matrix NEXUS
---------------------

Used for pairwise distance matrices. Blocks: **TAXA** (taxon labels) and
**Distances** (matrix).

**Classes using this subtype:** :class:`phylozoo.core.distance.DistanceMatrix` (default format).

Structure:

.. code-block:: text

   #NEXUS

   BEGIN Taxa;
       DIMENSIONS ntax=3;
       TAXLABELS
           A
           B
           C
       ;
   END;

   BEGIN Distances;
       DIMENSIONS ntax=3;
       FORMAT triangle=LOWER diagonal LABELS;
       MATRIX
       A 0.000000
       B 1.000000 0.000000
       C 2.000000 1.000000 0.000000
       ;
   END;

Example:

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np

   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   dm.save("distances.nexus")
   dm2 = DistanceMatrix.load("distances.nexus")
   nexus_str = dm.to_string(format="nexus", triangle="UPPER")

MSA NEXUS (Characters)
----------------------

Used for multiple sequence alignments. Blocks: **TAXA** and **CHARACTERS**.

**Classes using this subtype:** :class:`phylozoo.core.sequence.MSA`.

Structure:

.. code-block:: text

   #NEXUS

   BEGIN TAXA;
       DIMENSIONS ntax=2;
       TAXLABELS
           taxon1
           taxon2
       ;
   END;

   BEGIN CHARACTERS;
       DIMENSIONS nchar=8;
       FORMAT datatype=DNA missing=N gap=-;
       MATRIX
       taxon1    ACGTACGT
       taxon2    TGCAACGT
       ;
   END;

Example:

.. code-block:: python

   from phylozoo import MSA
   sequences = {"taxon1": "ACGTACGT", "taxon2": "TGCAACGT"}
   msa = MSA(sequences)
   msa.save("alignment.nexus", format="nexus")
   msa2 = MSA.load("alignment.nexus")

Split System NEXUS
------------------

Used for split systems (and weighted split systems). Blocks: **TAXA** and **SPLITS**.

**Classes using this subtype:** :class:`phylozoo.core.split.SplitSystem`,
:class:`phylozoo.core.split.WeightedSplitSystem` (default format).

Structure:

.. code-block:: text

   #NEXUS

   BEGIN TAXA;
       DIMENSIONS ntax=4;
       TAXLABELS
           A
           B
           C
           D
       ;
   END;

   BEGIN SPLITS;
       DIMENSIONS ntax=4;
       FORMAT labels=yes weights=yes;
       MATRIX
           1.0  A B | C D
           0.8  A C | B D
       ;
   END;

Example:

.. code-block:: python

   from phylozoo.core.split import SplitSystem, Split
   split1 = Split({'A', 'B'}, {'C', 'D'})
   split2 = Split({'A', 'C'}, {'B', 'D'})
   splits = SplitSystem([split1, split2])
   splits.save("splits.nexus")
   splits2 = SplitSystem.load("splits.nexus")

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`phylip` — Alternative distance matrix format
   * :doc:`fasta` — Alternative sequence format
   * :doc:`../../../core/distance` — Distance matrices
   * :doc:`../../../core/sequences` — MSA
   * :doc:`../../../core/splits/overview` — Split systems
