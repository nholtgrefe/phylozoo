Multiple Sequence Alignments
============================

The :mod:`phylozoo.core.sequence` module provides immutable, array-backed containers for multiple sequence
alignments (MSAs), along with comprehensive tools for distance computation, bootstrap
resampling, and other operations essential for phylogenetic analysis. MSAs serve as
the standard input format for many downstream algorithms, including quartet profile
generation and network reconstruction.

All classes and functions on this page can be imported from the core sequence module:

.. code-block:: python

   from phylozoo.core.sequence import *
   # or directly
   from phylozoo.core.sequence import MSA

Working with Multiple Sequence Alignments
------------------------------------------

The :class:`phylozoo.core.sequence.MSA` class is the canonical container for aligned
sequences in PhyloZoo. It provides an immutable, labeled, and read-only representation
that ensures data integrity throughout your analysis pipeline.

.. note::
   :class: dropdown

   **Implementation details**

   Multiple sequence alignments are designed for immutability and performance:

   - Sequences are stored internally as a NumPy ``int8`` array of shape ``(num_taxa, sequence_length)``
   - The default encoding maps nucleotides A/C/G/T to 0/1/2/3, with gaps/unknown characters as -1
   - Taxa are stored in canonical sorted order and exposed via immutable properties
   - To modify sequences, create a new ``MSA`` instance with the updated data
   - The alignment is validated for equal sequence lengths at construction time

   For implementation details, see :mod:`src/phylozoo/core/sequence/base.py`.

Creating an MSA
^^^^^^^^^^^^^^^^

MSAs can be created from dictionaries mapping taxon names to sequence strings:

.. code-block:: python

   from phylozoo.core.sequence import MSA

   # Create from sequence dictionary
   sequences = {
       "taxon1": "ACGTACGT",
       "taxon2": "ACGTACGT",
       "taxon3": "ACGTTTAA"
   }

   msa = MSA(sequences)

All sequences must have the same length. The constructor validates these properties
and raises :class:`phylozoo.utils.exceptions.PhyloZooValueError` if the input is invalid.

For performance-critical applications, you can also create MSAs directly from pre-encoded arrays:

.. code-block:: python

   import numpy as np

   # Create from coded array (efficient for internal operations)
   coded_array = np.array([[0, 1, 2, 3, 0, 1, 2, 3],  # taxon1
                          [0, 1, 2, 3, 0, 1, 2, 3],  # taxon2
                          [0, 1, 2, 3, 3, 3, 0, 0]], # taxon3
                         dtype=np.int8)
   taxa_order = ("taxon1", "taxon2", "taxon3")

   msa = MSA.from_coded_array(coded_array, taxa_order)

Accessing Sequences
^^^^^^^^^^^^^^^^^^^

MSAs provide several methods for accessing sequences and metadata:

.. code-block:: python

   # Get sequence for a specific taxon
   sequence = msa.get_sequence("taxon1")  # Returns: "ACGTACGT"

   # Check if taxon exists
   exists = "taxon1" in msa  # Returns: True

   # Access metadata
   taxa = msa.taxa          # Returns: frozenset of taxon names
   taxa_order = msa.taxa_order  # Returns: tuple with canonical order
   length = msa.sequence_length # Returns: 8
   num_taxa = msa.num_taxa      # Returns: 3

   # Access internal representation
   coded = msa.coded_array  # Returns: read-only numpy array

These methods provide safe access to sequences without exposing mutable internals,
maintaining immutability guarantees.

File Input/Output
^^^^^^^^^^^^^^^^^

MSAs support reading and writing in multiple phylogenetic formats:

- **FASTA** (default): Standard sequence alignment format
- **NEXUS**: Comprehensive phylogenetic data format

.. code-block:: python

   # Load from file (auto-detects format by extension)
   msa = MSA.load("alignment.fasta")

   # Load with explicit format
   msa = MSA.load("alignment.nexus", format="nexus")

   # Save to file
   msa.save("output.fasta")

.. seealso::
   The I/O system uses the :class:`phylozoo.utils.io.IOMixin` interface, providing
   consistent file handling across PhyloZoo classes. For details on the I/O system,
   see the :doc:`I/O documentation <../utils/io>`. For specific information about
   supported file formats and parameter options for MSAs, see the
   :mod:`API reference <phylozoo.core.sequence.io>`.

Bootstrap Resampling
---------------------

The sequence module provides functions for bootstrap resampling, which is essential
for assessing the statistical support of phylogenetic inferences.

Basic Bootstrap
^^^^^^^^^^^^^^^

The :func:`phylozoo.core.sequence.bootstrap` function generates bootstrap replicates
by sampling alignment columns with replacement:

.. code-block:: python

   from phylozoo.core.sequence import bootstrap

   # Generate 100 bootstrap replicates
   for replicate in bootstrap(msa, n_bootstrap=100, seed=42):
       # Each replicate is a new MSA with resampled columns
       print(f"Replicate has {replicate.sequence_length} columns")

The ``length`` parameter controls how many columns to sample (defaults to full alignment
length), and ``seed`` ensures reproducible results for testing and debugging.

Gene-Based Bootstrap
^^^^^^^^^^^^^^^^^^^^

For multi-gene alignments, the :func:`phylozoo.core.sequence.bootstrap_per_gene`
function resamples columns within each gene separately:

.. code-block:: python

   from phylozoo.core.sequence import bootstrap_per_gene

   # Define gene boundaries (lengths must sum to total alignment length)
   gene_lengths = [100, 200, 150]  # Three genes

   for replicate in bootstrap_per_gene(msa, gene_lengths, n_bootstrap=100, seed=42):
       # Columns are resampled within each gene independently
       print(f"Gene-based replicate: {replicate.sequence_length} columns")

Distance Computation
--------------------

The sequence module provides efficient functions for computing evolutionary distances
from multiple sequence alignments.

Hamming Distance
^^^^^^^^^^^^^^^^

The :func:`phylozoo.core.sequence.hamming_distances` function computes normalized
Hamming distances between all pairs of sequences:

.. math::

   d(i,j) = \frac{1}{L} \sum_{k=1}^{L} \mathbf{1}_{s_i[k] \neq s_j[k]}

where :math:`L` is the number of valid (non-gap, non-unknown) positions, and the
indicator function equals 1 when sequences differ at position :math:`k`.

.. code-block:: python

   from phylozoo.core.sequence import hamming_distances
   from phylozoo.core.distance.classifications import is_metric

   # Compute distance matrix
   distance_matrix = hamming_distances(msa)

   # Check if result is a proper metric
   if is_metric(distance_matrix):
       print("Hamming distances form a metric")

The function excludes positions where either sequence contains gaps (-) or unknown
characters (N), focusing only on positions where both sequences have valid nucleotides.
The implementation uses vectorized NumPy operations for efficient computation on
large alignments.

See Also
--------

- :doc:`API Reference <../../api/core/sequence>` - Complete function signatures and detailed examples
- :mod:`phylozoo.core.sequence.io` - MSA I/O format details and parameter options
- :doc:`Distance Matrices <distance>` - Working with distance matrices computed from alignments