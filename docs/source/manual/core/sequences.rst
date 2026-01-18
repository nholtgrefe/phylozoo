Sequences
=========

This page covers working with multiple sequence alignments (MSAs) in PhyloZoo, including 
basic operations, bootstrapping, and advanced array operations.

**I/O Formats**: FASTA (default, extensions: ``.fasta``, ``.fa``, ``.fas``), NEXUS (extensions: ``.nexus``, ``.nex``, ``.nxs``). 
See :doc:`I/O <../io>` for details.

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

See the :doc:`I/O <../reference/io>` page for supported formats and detailed I/O operations.

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

Per-Gene Bootstrapping
-----------------------

When working with multi-gene alignments, you may want to bootstrap within gene 
boundaries:

.. code-block:: python

   from phylozoo.core.sequence import bootstrap_per_gene
   
   # Define gene boundaries
   gene_lengths = [100, 150, 200]  # Lengths of each gene
   
   # Bootstrap per gene
   for bootstrapped_msa in bootstrap_per_gene(
       msa, 
       gene_lengths, 
       n_bootstrap=100, 
       seed=42
   ):
       # Process each bootstrap replicate
       # Columns are resampled within each gene separately
       pass

This ensures that bootstrap replicates maintain gene structure, which is important 
for multi-gene phylogenetic analyses.

Efficient Array Operations
--------------------------

MSAs store sequences internally as NumPy arrays for efficient computation. You can 
work directly with the coded array representation:

.. code-block:: python

   # Access the internal NumPy array
   coded_array = msa.coded_array  # Shape: (num_taxa, sequence_length), dtype: int8
   
   # Create MSA from coded array (efficient, avoids encoding/decoding)
   from phylozoo import MSA
   new_msa = MSA.from_coded_array(coded_array, msa.taxa_order)

This is particularly useful for custom operations that work directly with the 
numerical representation.

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

**Advanced Functions:**

* **bootstrap_per_gene(msa, gene_lengths, n_bootstrap, seed=None)** - Generate 
  bootstrap replicates with per-gene resampling. Columns are resampled separately 
  within each gene boundary. Gene lengths must sum to the total sequence length. 
  Returns iterator of MSA objects. See :func:`phylozoo.core.sequence.bootstrap_per_gene` 
  for details.

**Basic MSA Methods:**

* **num_taxa** - Number of taxa in the alignment
* **sequence_length** - Length of sequences (all sequences have same length)
* **taxa** - Frozen set of all taxon names
* **taxa_order** - Tuple of taxon names in canonical order
* **get_sequence(taxon)** - Get sequence string for a taxon
* **save(filename)** - Save MSA to file
* **load(filename)** - Load MSA from file (class method)

**Advanced MSA Methods:**

* **coded_array** - Access the internal NumPy array representation. Shape is 
  (num_taxa, sequence_length) with dtype int8. Each element is an encoded character 
  value. Read-only property.

* **from_coded_array(coded_array, taxa_order)** - Create MSA directly from coded 
  NumPy array. Efficient method that avoids encoding/decoding overhead. Useful for 
  operations that work directly with arrays. Class method. See 
  :meth:`phylozoo.core.sequence.MSA.from_coded_array` for details.

.. note::
   MSAs are immutable. To modify an MSA, create a new instance with the desired 
   changes.

.. note::
   The coded array uses DEFAULT_NUCLEOTIDE_ENCODING: A/C/G/T = 0/1/2/3, gaps/unknown = -1.

.. tip::
   Use ``from_coded_array`` when you're doing many bootstrap operations or custom 
   array manipulations. It's much faster than creating MSAs from strings.
