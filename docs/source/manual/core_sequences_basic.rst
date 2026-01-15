Sequences (Basic)
=================

This page covers basic operations for working with multiple sequence alignments (MSAs). 
For advanced sequence operations like per-gene bootstrapping, see 
:doc:`Sequences (Advanced) <core_sequences_advanced>`.

Creating MSAs
-------------

Create MSAs from dictionaries of sequences:

.. code-block:: python

   from phylozoo import MSA
   
   sequences = {
       "taxon1": "ACGTACGT",
       "taxon2": "ACGTACGT",
       "taxon3": "ACGTTTAA"
   }
   
   msa = MSA(sequences)
   print(f"Number of taxa: {msa.num_taxa}")
   print(f"Sequence length: {msa.sequence_length}")

All sequences must have the same length. The MSA class validates this automatically.

Loading and Saving
------------------

MSAs can be loaded from and saved to files:

.. code-block:: python

   # Save MSA
   msa.save("alignment.fasta")
   
   # Load MSA
   msa = MSA.load("alignment.fasta")
   msa = MSA.load("alignment.nexus", format="nexus")

See the :doc:`I/O <io>` page for supported formats and detailed I/O operations.

Accessing Sequences
-------------------

Access sequences by taxon name:

.. code-block:: python

   # Get sequence for a taxon
   sequence = msa.get_sequence("taxon1")
   
   # Check if taxon exists
   if "taxon1" in msa:
       print("Taxon found")
   
   # Iterate over all taxa
   for taxon in msa.taxa:
       seq = msa.get_sequence(taxon)
       print(f"{taxon}: {seq}")

Basic Bootstrapping
-------------------

Generate bootstrap replicates:

.. code-block:: python

   from phylozoo.core.sequence import bootstrap
   
   # Generate 100 bootstrap replicates
   for bootstrapped_msa in bootstrap(msa, n_bootstrap=100, seed=42):
       # Process each bootstrap replicate
       distances = hamming_distances(bootstrapped_msa)
       # ... use distances for analysis

Distance Calculations
----------------------

Compute distance matrices from MSAs:

.. code-block:: python

   from phylozoo.core.sequence import hamming_distances
   
   # Compute normalized Hamming distances
   distance_matrix = hamming_distances(msa)
   
   # Access distances
   dist = distance_matrix.get_distance("taxon1", "taxon2")

The ``hamming_distances`` function computes normalized Hamming distances, excluding 
positions with gaps or unknown characters.

Available Functions
-------------------

**Classes:**

* **MSA** - Multiple Sequence Alignment class. Immutable class for working with 
  aligned sequences. Sequences are stored efficiently as NumPy arrays. Supports 
  bootstrapping and distance calculations. See :class:`phylozoo.core.sequence.MSA` 
  for full API.

**Basic Functions:**

* **bootstrap(msa, n_bootstrap, length=None, seed=None)** - Generate bootstrap 
  replicates by resampling columns with replacement. Returns iterator of MSA objects. 
  The length parameter allows resampling a specific number of columns. See 
  :func:`phylozoo.core.sequence.bootstrap` for details.

* **hamming_distances(msa)** - Compute normalized Hamming distances between all pairs 
  of taxa. Excludes positions with gaps or unknown characters. Uses vectorized NumPy 
  operations for efficiency. Returns DistanceMatrix. See 
  :func:`phylozoo.core.sequence.hamming_distances` for details.

**Basic MSA Methods:**

* **num_taxa** - Number of taxa in the alignment
* **sequence_length** - Length of sequences (all sequences have same length)
* **taxa** - Frozen set of all taxon names
* **taxa_order** - Tuple of taxon names in canonical order
* **get_sequence(taxon)** - Get sequence string for a taxon
* **save(filename)** - Save MSA to file
* **load(filename)** - Load MSA from file (class method)

.. note::
   MSAs are immutable. To modify an MSA, create a new instance with the desired 
   changes.

.. tip::
   Use bootstrapping to assess confidence in your analyses. Generate multiple 
   bootstrap replicates and analyze each one to see how robust your results are.

.. seealso::
   For per-gene bootstrapping and advanced sequence operations, see 
   :doc:`Sequences (Advanced) <core_sequences_advanced>`.
