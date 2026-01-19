Sequence Alignments
===================

The sequence module provides immutable, array-backed containers for multiple sequence
alignments (MSAs), together with utilities for distance computation and bootstrap
resampling commonly used in phylogenetic analyses.

All classes and functions on this page can be imported from the core sequence module:

.. code-block:: python

   from phylozoo.core import sequence
   # or directly
   from phylozoo.core.sequence import MSA

Working with MSAs
-----------------

The :class:`phylozoo.core.sequence.MSA` class is the canonical, immutable container for
aligned sequences in PhyloZoo. It stores sequences in a canonical order and uses an
internal coded NumPy array representation for efficient processing.

.. note::
   :class: dropdown

   **Implementation details**

   - Sequences are stored internally as a NumPy ``int8`` array (`_coded_array`) of shape
     ``(num_taxa, sequence_length)``.
   - The default nucleotide encoding maps A/C/G/T to 0/1/2/3, with gaps/unknown as -1.
   - Taxa are stored in canonical sorted order and exposed via ``taxa_order`` and
     ``taxa`` properties.
   - I/O uses the :class:`~phylozoo.utils.io.IOMixin` common interface; supported
     formats include FASTA (default) and NEXUS. See
     :mod:`src/phylozoo/core/sequence/base.py` for details.

Creating an MSA
^^^^^^^^^^^^^^^

Create a :class:`phylozoo.core.sequence.MSA` from a dictionary mapping taxon names to
sequence strings. All sequences must have the same length.

.. code-block:: python

   from phylozoo.core.sequence import MSA

   sequences = {
       "taxon1": "ACGTACGT",
       "taxon2": "ACGTACGT",
       "taxon3": "ACGTTTAA"
   }

   msa = MSA(sequences)  # taxa_order becomes a tuple

The constructor validates sequence lengths and types, and will raise
:class:`phylozoo.utils.exceptions.PhyloZooValueError` or
:class:`phylozoo.utils.exceptions.PhyloZooTypeError` for invalid inputs.

File Input/Output
^^^^^^^^^^^^^^^^^

Use the I/O helpers to read and write MSAs in common formats.

.. code-block:: python

   # Load (auto-detect by extension)
   msa = MSA.load("alignment.fasta")

   # Load with explicit format
   msa = MSA.load("alignment.nexus", format="nexus")

   # Save
   msa.save("out.fasta")

The I/O system delegates to the project's format registry and readers via
:class:`~phylozoo.utils.io.IOMixin`. For details about supported formats and
additional parameters see the global I/O docs.

Class methods / Accessing sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Common accessors and conveniences:

.. code-block:: python

   msa.get_sequence("taxon1")   # -> str or None
   "taxon1" in msa              # -> bool
   msa.taxa                      # -> frozenset of taxon names
   msa.taxa_order                # -> tuple with canonical order
   msa.sequence_length           # -> int
   msa.num_taxa                  # -> int
   msa.coded_array               # -> numpy array (num_taxa, sequence_length)

Important class methods:

- ``MSA.from_coded_array(coded_array, taxa_order)`` — Create an MSA directly from a
  coded NumPy array (dtype ``int8``). Efficient for operations like bootstrapping and
  avoids repeated encoding/decoding.

Operations
----------

The :mod:`phylozoo.core.sequence` module provides utilities for bootstrap resampling
and distance computation from MSAs. These operations are commonly used to prepare
inputs for downstream inference algorithms.

Bootstrap resampling
^^^^^^^^^^^^^^^^^^^^

- ``bootstrap(msa, n_bootstrap, length=None, seed=None)`` — Generate bootstrap
  replicates by sampling columns with replacement. ``length`` sets the number of
  columns sampled (defaults to alignment length), and ``seed`` ensures reproducibility.

- ``bootstrap_per_gene(msa, gene_lengths, n_bootstrap, seed=None)`` — For multi-gene
  alignments, resample columns within each gene separately (``gene_lengths`` must sum
  to the total sequence length).

Example:

.. code-block:: python

   from phylozoo.core.sequence import bootstrap, hamming_distances

   for rep in bootstrap(msa, n_bootstrap=100, seed=42):
       dm = hamming_distances(rep)
       # use dm in analysis

Distance computation
^^^^^^^^^^^^^^^^^^^^

- ``hamming_distances(msa)`` — Compute normalized Hamming distances between all pairs
  of taxa. Positions where either sequence contains a gap/unknown (encoded as ``-1``)
  are excluded from pairwise comparisons. Returns a
  :class:`phylozoo.core.distance.DistanceMatrix` suitable for classification and
  inference.

  The normalization excludes positions where either sequence has a gap (-) or unknown
  character (N). Only positions where both sequences have valid nucleotides (A, C, G, T)
  are considered. The implementation uses vectorized NumPy operations for efficiency.

Performance and testing notes
----------------------------

- ``MSA.coded_array`` is the canonical internal representation. For performance-sensitive
  code, operate directly on the coded array and then re-create MSAs using
  ``MSA.from_coded_array`` to avoid repeated encoding/decoding overhead.
- When adding algorithms that manipulate the coded array, benchmark and document
  speedups (use ``sandbox/findings.txt``). Prefer vectorized NumPy or Numba-accelerated
  implementations for heavy computations.
- For reproducible tests involving randomness (bootstrapping), pin the seed and assert
  deterministic properties (e.g., that replicates are identical for the same seed).
- Record micro-benchmarks and experimental findings in ``sandbox/findings.txt``.

Notes and tips
-------------

- MSAs are immutable by design; to change values create a new :class:`MSA` instance
  to avoid shared-state bugs.
- Prefer ``from_coded_array`` for heavy array workloads to avoid repeated encoding.
- Validate MSAs early (sequence lengths, valid characters) to produce clearer error
  messages and stable downstream behavior.
- When adding examples to docstrings, mirror them in tests under ``tests/`` as per
  project conventions.
