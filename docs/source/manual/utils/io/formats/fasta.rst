FASTA Format
============

FASTA is a simple text-based format for sequence data, widely used in bioinformatics.
Each sequence has a header line starting with ``>`` followed by one or more lines of
sequence characters.

**Classes using FASTA format:** :class:`phylozoo.core.sequence.MSA` (default format).

**File extensions:** ``.fasta``, ``.fa``, ``.fas``

Format Structure
----------------

.. code-block:: text

   >taxon1
   ACGTACGT
   >taxon2
   TGCAACGT
   >taxon3
   AAAAACGT

Example
-------

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

.. seealso::
   * :doc:`../operations` — Save/load and format detection
   * :doc:`nexus` — NEXUS Characters block for alignments
   * :doc:`../../../core/sequences` — MSA and sequences
