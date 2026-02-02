sequence
========

The sequence module provides classes and functions for working with biological sequences
and multiple sequence alignments (MSAs). Sequences are stored internally as NumPy arrays
for efficient computation.

Main Classes
------------

.. automodule:: phylozoo.core.sequence.base
   :members:
   :show-inheritance:

Distance Computation
--------------------

Compute distance matrices from sequence alignments using normalized Hamming distances.
The Hamming distance between two sequences is the number of positions where they differ,
normalized by the number of valid positions (excluding gaps and unknown characters).

.. automodule:: phylozoo.core.sequence.distances
   :members:
   :show-inheritance:

Bootstrap Functions
------------------

Generate bootstrap replicates of sequence alignments for statistical analysis.

.. automodule:: phylozoo.core.sequence.bootstrap
   :members:
   :show-inheritance:

I/O Support
-----------

MSAs support reading and writing in FASTA and NEXUS formats.

.. automodule:: phylozoo.core.sequence.io
   :members:
   :show-inheritance:
