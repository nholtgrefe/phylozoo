Sequences (Advanced)
====================

This page covers advanced sequence operations including per-gene bootstrapping and 
efficient array operations. For basic MSA operations, see 
:doc:`Sequences (Basic) <core_sequences_basic>`.

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
---------------------------

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

**Advanced Functions:**

* **bootstrap_per_gene(msa, gene_lengths, n_bootstrap, seed=None)** - Generate 
  bootstrap replicates with per-gene resampling. Columns are resampled separately 
  within each gene boundary. Gene lengths must sum to the total sequence length. 
  Returns iterator of MSA objects. See :func:`phylozoo.core.sequence.bootstrap_per_gene` 
  for details.

**Advanced MSA Methods:**

* **coded_array** - Access the internal NumPy array representation. Shape is 
  (num_taxa, sequence_length) with dtype int8. Each element is an encoded character 
  value. Read-only property.

* **from_coded_array(coded_array, taxa_order)** - Create MSA directly from coded 
  NumPy array. Efficient method that avoids encoding/decoding overhead. Useful for 
  operations that work directly with arrays. Class method. See 
  :meth:`phylozoo.core.sequence.MSA.from_coded_array` for details.

.. note::
   The coded array uses DEFAULT_NUCLEOTIDE_ENCODING: A/C/G/T = 0/1/2/3, gaps/unknown = -1.

.. tip::
   Use ``from_coded_array`` when you're doing many bootstrap operations or custom 
   array manipulations. It's much faster than creating MSAs from strings.

.. warning::
   When using ``from_coded_array``, ensure the coded array uses the correct encoding. 
   The method validates the array but assumes DEFAULT_NUCLEOTIDE_ENCODING.
