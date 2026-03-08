FASTA
=====

FASTA is a simple text-based format for sequence data, widely used in bioinformatics.
Each sequence has a header line starting with ``>`` followed by one or more lines of
sequence characters.

.. seealso::
   `FASTA format <https://en.wikipedia.org/wiki/FASTA_format>`_ — Wikipedia

Classes and extensions
----------------------

**Classes:** :class:`~phylozoo.core.sequence.base.MSA` (default format)

**File extensions:** ``.fasta``, ``.fa``, ``.fas``

Structure
---------

.. code-block:: text

   >taxon1
   ACGTACGT
   >taxon2
   TGCAACGT
   >taxon3
   AAAAACGT

Examples
--------

.. code-block:: python

   from phylozoo import MSA
   sequences = {
       "taxon1": "ACGTACGT",
       "taxon2": "TGCAACGT",
       "taxon3": "AAAAACGT"
   }
   msa = MSA(sequences)
   msa.save("alignment.fasta")
   msa2 = MSA.load("alignment.fasta")
   fasta_str = msa.to_string(format="fasta", line_length=60)

See also
--------

- :doc:`../operations` — Save/load and format detection
- :doc:`nexus` — NEXUS Characters block for alignments
- :doc:`../../../core/sequences` — MSA and sequences
