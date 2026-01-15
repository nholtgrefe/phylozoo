Input/Output (I/O)
==================

PhyloZoo provides a unified I/O system for reading and writing phylogenetic data in 
various formats. All core classes inherit from ``IOMixin``, which provides consistent 
methods for file operations.

Basic I/O Operations
--------------------

All I/O-capable classes support these methods:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork, MSA, DistanceMatrix
   
   # Save to file (auto-detects format from extension)
   network.save("network.enewick")
   msa.save("alignment.fasta")
   dm.save("distances.nexus")
   
   # Load from file (auto-detects format)
   network = DirectedPhyNetwork.load("network.enewick")
   msa = MSA.load("alignment.fasta")
   dm = DistanceMatrix.load("distances.nexus")
   
   # Convert to string
   enewick_string = network.to_string(format="enewick")
   
   # Parse from string
   network = DirectedPhyNetwork.from_string(enewick_string, format="enewick")
   
   # Convert between formats
   MSA.convert("input.fasta", "output.nexus")

Supported Formats
-----------------

Network Formats
~~~~~~~~~~~~~~~

**eNewick** (Extended Newick)
   Default format for networks. Supports both trees and networks with hybrid nodes. 
   See :cite:`Cardona2008` for the format specification.
   
   Extensions: ``.enewick``, ``.enewick``
   
   Example:
   
   .. code-block:: python
      
      network = DirectedPhyNetwork.load("network.enewick")
      network.save("output.enewick")

**DOT** (Graphviz)
   Graphviz format for networks. Useful for visualization and integration with 
   Graphviz tools.
   
   Extensions: ``.dot``, ``.gv``
   
   Example:
   
   .. code-block:: python
      
      network.save("network.dot", format="dot")

**Edge List**
   Simple edge list format for basic network structures.
   
   Extensions: ``.edgelist``, ``.edges``

Sequence Formats
~~~~~~~~~~~~~~~~

**FASTA**
   Standard FASTA format for sequence alignments. Default format for MSAs.
   
   Extensions: ``.fasta``, ``.fa``, ``.fas``
   
   Example:
   
   .. code-block:: python
      
      msa = MSA.load("alignment.fasta")
      msa.save("output.fasta")

**NEXUS**
   NEXUS format for sequence alignments.
   
   Extensions: ``.nexus``, ``.nex``, ``.nxs``
   
   Example:
   
   .. code-block:: python
      
      msa = MSA.load("alignment.nexus", format="nexus")

Distance Matrix Formats
~~~~~~~~~~~~~~~~~~~~~~~

**NEXUS**
   NEXUS format for distance matrices.
   
   Extensions: ``.nexus``, ``.nex``, ``.nxs``
   
   Example:
   
   .. code-block:: python
      
      dm = DistanceMatrix.load("distances.nexus")

**PHYLIP**
   PHYLIP format for distance matrices. Supports both lower and upper triangle formats.
   
   Extensions: ``.phy``, ``.phylip``
   
   Example:
   
   .. code-block:: python
      
      dm = DistanceMatrix.load("distances.phy", format="phylip")

**CSV**
   Comma-separated values format for distance matrices.
   
   Extensions: ``.csv``
   
   Example:
   
   .. code-block:: python
      
      dm = DistanceMatrix.load("distances.csv", format="csv")

Split System Formats
~~~~~~~~~~~~~~~~~~~~

**NEXUS**
   NEXUS format for split systems.
   
   Extensions: ``.nexus``, ``.nex``, ``.nxs``

Format Detection
----------------

PhyloZoo automatically detects file formats based on file extensions. If the extension 
is ambiguous or you want to specify the format explicitly, use the ``format`` parameter:

.. code-block:: python

   # Explicit format specification
   network = DirectedPhyNetwork.load("network.txt", format="enewick")
   msa = MSA.load("data.txt", format="fasta")

Error Handling
--------------

I/O operations raise appropriate exceptions when files cannot be read or written:

* ``PhyloZooIOError``: General I/O errors
* ``PhyloZooParseError``: Parsing errors (malformed files)
* ``FileNotFoundError``: File not found errors

Example:

.. code-block:: python

   from phylozoo import PhyloZooParseError
   
   try:
       network = DirectedPhyNetwork.load("malformed.enewick")
   except PhyloZooParseError as e:
       print(f"Parse error: {e}")
