NEXUS
=====

The NEXUS format is a flexible, block-structured format used widely in phylogenetics.
A NEXUS file starts with ``#NEXUS`` and contains one or more blocks, each of the form
``BEGIN blockname; ... END;``. PhyloZoo uses different block types for different
data kinds (distance matrices, sequence alignments, split systems). The same file
extensions are used for all NEXUS subtypes; the class you load or save with
determines which block structure is expected.

.. seealso::
   `Nexus file <https://en.wikipedia.org/wiki/Nexus_file>`_ — Wikipedia

Classes and extensions
----------------------

**File extensions:** ``.nexus``, ``.nex``, ``.nxs``

**Classes (by subtype):** :class:`~phylozoo.core.distance.base.DistanceMatrix` (DISTANCES),
:class:`~phylozoo.core.sequence.base.MSA` (CHARACTERS),
:class:`~phylozoo.core.split.splitsystem.SplitSystem` and
:class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem` (SPLITS).

Structure
---------

All NEXUS files begin with the ``#NEXUS`` token. Data are organized in blocks.
Each block has ``BEGIN blockname;``, block-specific commands and a ``MATRIX`` or
data section, and ``END;``. The **TAXA** block (with ``TAXLABELS``) is shared across
subtypes to define taxon labels. The data block (DISTANCES, CHARACTERS, or SPLITS)
holds the actual data.

.. code-block:: text

   #NEXUS

   BEGIN TAXA;
       DIMENSIONS ntax=N;
       TAXLABELS
           taxon1
           taxon2
       ;
   END;

   BEGIN SomeBlock;
       ...
   END;

DISTANCES
^^^^^^^^^

Blocks: **TAXA** and **DISTANCES**. The DISTANCES block contains a lower or upper
triangular matrix with optional ``FORMAT triangle=...``.

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

   BEGIN DISTANCES;
       DIMENSIONS ntax=3;
       FORMAT triangle=LOWER diagonal LABELS;
       MATRIX
       A 0.000000
       B 1.000000 0.000000
       C 2.000000 1.000000 0.000000
       ;
   END;

CHARACTERS
^^^^^^^^^^^

Blocks: **TAXA** and **CHARACTERS**. The Characters block has ``DIMENSIONS nchar=...``,
``FORMAT datatype=...``, and a ``MATRIX`` with aligned sequences.

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

SPLITS
^^^^^^

Blocks: **TAXA** and **SPLITS**. The SPLITS block has ``FORMAT labels=yes weights=yes``
and a ``MATRIX`` with splits in ``A B | C D`` notation, optionally with weights.

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

Examples
--------

**Distance matrix:**

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np

   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   dm.save("distances.nexus")
   dm2 = DistanceMatrix.load("distances.nexus")
   nexus_str = dm.to_string(format="nexus", triangle="UPPER")

**MSA:**

.. code-block:: python

   from phylozoo import MSA
   sequences = {"taxon1": "ACGTACGT", "taxon2": "TGCAACGT"}
   msa = MSA(sequences)
   msa.save("alignment.nexus", format="nexus")
   msa2 = MSA.load("alignment.nexus")

**Split system:**

.. code-block:: python

   from phylozoo.core.split import SplitSystem, Split
   split1 = Split({'A', 'B'}, {'C', 'D'})
   split2 = Split({'A', 'C'}, {'B', 'D'})
   splits = SplitSystem([split1, split2])
   splits.save("splits.nexus")
   splits2 = SplitSystem.load("splits.nexus")

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`phylip` — Alternative distance matrix format
- :doc:`fasta` — Alternative sequence format
- :doc:`../../../core/distance` — Distance matrices
- :doc:`../../../core/sequences` — MSA
- :doc:`../../../core/splits/overview` — Split systems
