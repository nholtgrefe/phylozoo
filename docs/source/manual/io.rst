Input/Output (I/O)
==================

PhyloZoo provides a unified I/O system for reading and writing phylogenetic data in 
various formats. All core classes inherit from ``IOMixin``, which provides consistent 
methods for file operations.

I/O System Overview
-------------------

All I/O-capable classes in PhyloZoo use the ``IOMixin`` protocol, which provides a 
consistent interface for loading, saving, and converting between formats. The system 
automatically detects file formats based on extensions and supports explicit format 
specification.

Basic I/O Operations
-------------------

All I/O-capable classes support these methods:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork, MSA, DistanceMatrix
   
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

Format Detection
----------------

PhyloZoo automatically detects file formats based on file extensions. If the extension 
is ambiguous or you want to specify the format explicitly, use the ``format`` parameter:

.. code-block:: python

   # Explicit format specification
   network = DirectedPhyNetwork.load("network.txt", format="enewick")
   msa = MSA.load("data.txt", format="fasta")

Supported Formats by Class
--------------------------

**DirectedPhyNetwork**
   * **eNewick** (default): Extended Newick format for networks. Extensions: ``.enewick``, ``.eNewick``, ``.enwk``, ``.nwk``, ``.newick``
   * **DOT**: Graphviz format. Extensions: ``.dot``, ``.gv``

**SemiDirectedPhyNetwork**
   * **Newick** (default): Standard Newick format. Extensions: ``.nwk``, ``.newick``, ``.enewick``, ``.eNewick``, ``.enw``
   * **PhyloZoo-DOT**: Custom DOT format for semi-directed networks. Extension: ``.pzdot``

**MSA**
   * **FASTA** (default): Standard FASTA format. Extensions: ``.fasta``, ``.fa``, ``.fas``
   * **NEXUS**: NEXUS format for sequence alignments. Extensions: ``.nexus``, ``.nex``, ``.nxs``

**DistanceMatrix**
   * **NEXUS** (default): NEXUS format for distance matrices. Extensions: ``.nexus``, ``.nex``, ``.nxs``
   * **PHYLIP**: PHYLIP format. Extensions: ``.phy``, ``.phylip``
   * **CSV**: Comma-separated values format. Extension: ``.csv``

**SplitSystem** and **WeightedSplitSystem**
   * **NEXUS** (default): NEXUS format for split systems. Extensions: ``.nexus``, ``.nex``, ``.nxs``

**SqQuartetProfileSet**
   * **PZ** (default): Custom JSON-based format. Extension: ``.pz``

Format Details
--------------

Network Formats
~~~~~~~~~~~~~~~~

**eNewick** (Extended Newick)
   Default format for ``DirectedPhyNetwork``. Supports both trees and networks with hybrid nodes. 
   See :cite:`Cardona2008` for the format specification. Preserves branch lengths, bootstrap 
   support, and gamma values for hybrid edges.
   
   Example:
   
   .. code-block:: python
      
      network = DirectedPhyNetwork.load("network.enewick")
      network.save("output.enewick")

**Newick**
   Standard Newick format for ``SemiDirectedPhyNetwork``. Supports trees and networks with 
   hybrid nodes. Undirected edges are represented in a way compatible with standard Newick parsers.
   
   Example:
   
   .. code-block:: python
      
      network = SemiDirectedPhyNetwork.load("network.newick")
      network.save("output.newick")

**DOT** (Graphviz)
   Graphviz format for networks. Useful for visualization and integration with Graphviz tools. 
   Preserves node labels and edge attributes.
   
   Example:
   
   .. code-block:: python
      
      network.save("network.dot", format="dot")

**PhyloZoo-DOT**
   Custom DOT format specifically designed for ``SemiDirectedPhyNetwork``. Preserves the 
   distinction between directed and undirected edges.
   
   Example:
   
   .. code-block:: python
      
      network.save("network.pzdot", format="phylozoo-dot")

Sequence Formats
~~~~~~~~~~~~~~~~

**FASTA**
   Standard FASTA format for sequence alignments. Default format for MSAs. Each sequence is 
   represented as a header line (starting with ``>``) followed by sequence data.
   
   Example:
   
   .. code-block:: python
      
      msa = MSA.load("alignment.fasta")
      msa.save("output.fasta")

**NEXUS**
   NEXUS format for sequence alignments. Supports metadata and multiple data blocks.
   
   Example:
   
   .. code-block:: python
      
      msa = MSA.load("alignment.nexus", format="nexus")

Distance Matrix Formats
~~~~~~~~~~~~~~~~~~~~~~~

**NEXUS**
   NEXUS format for distance matrices. Default format for ``DistanceMatrix``. Supports metadata 
   and multiple data blocks.
   
   Example:
   
   .. code-block:: python
      
      dm = DistanceMatrix.load("distances.nexus")

**PHYLIP**
   PHYLIP format for distance matrices. Supports both lower and upper triangle formats. 
   Compatible with standard PHYLIP tools.
   
   Example:
   
   .. code-block:: python
      
      dm = DistanceMatrix.load("distances.phy", format="phylip")

**CSV**
   Comma-separated values format for distance matrices. Easy to read and edit in spreadsheet 
   applications.
   
   Example:
   
   .. code-block:: python
      
      dm = DistanceMatrix.load("distances.csv", format="csv")

Split System Formats
~~~~~~~~~~~~~~~~~~~~

**NEXUS**
   NEXUS format for split systems. Supports both unweighted and weighted split systems. 
   Default format for ``SplitSystem`` and ``WeightedSplitSystem``.
   
   Example:
   
   .. code-block:: python
      
      splits = SplitSystem.load("splits.nexus")

Error Handling
--------------

I/O operations raise appropriate exceptions when files cannot be read or written:

* ``PhyloZooIOError``: General I/O errors
* ``PhyloZooParseError``: Parsing errors (malformed files)
* ``PhyloZooFormatError``: Format-related errors (unsupported format, missing format handler)
* ``FileNotFoundError``: File not found errors

Example:

.. code-block:: python

   from phylozoo import PhyloZooParseError
   
   try:
       network = DirectedPhyNetwork.load("malformed.enewick")
   except PhyloZooParseError as e:
       print(f"Parse error: {e}")

Best Practices
--------------

1. **Use appropriate file extensions**: PhyloZoo auto-detects formats from extensions. 
   Use standard extensions (e.g., ``.enewick``, ``.fasta``, ``.nexus``) for best results.

2. **Specify format explicitly when needed**: If the file extension is ambiguous or you're 
   working with non-standard extensions, specify the format explicitly.

3. **Handle exceptions**: Always use try/except blocks for I/O operations to handle 
   parsing errors and missing files gracefully.

4. **Preserve metadata**: Some formats preserve more metadata than others. Use formats 
   that preserve the attributes you need (e.g., eNewick for networks with branch lengths 
   and gamma values).

5. **Convert between formats**: Use the ``convert`` class method to convert between formats 
   without loading into memory first (useful for large files).

.. seealso::
   For more information on specific classes and their I/O capabilities, see:
   * :doc:`Networks (Basic) <core/networks/basic>` for network I/O
   * :doc:`Sequences <core/sequences>` for MSA I/O
   * :doc:`Distance Matrices <core/distance>` for distance matrix I/O
