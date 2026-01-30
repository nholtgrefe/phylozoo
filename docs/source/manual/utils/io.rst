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

Supported File Formats
======================

NEXUS Format
------------

The NEXUS format is a flexible, block-structured format used by multiple PhyloZoo classes.
It consists of blocks (e.g., ``BEGIN TAXA; ... END;``) that contain different types of data.

**Classes using NEXUS format:**

* **DistanceMatrix** (default format)
* **MSA** (Multiple Sequence Alignments)
* **SplitSystem** and **WeightedSplitSystem** (default format)

**File Extensions:** ``.nexus``, ``.nex``, ``.nxs``

**Distance Matrix NEXUS Format**

The NEXUS format for distance matrices includes a TAXA block and a DISTANCES block:

.. code-block:: text

   #NEXUS

   BEGIN Taxa;
       DIMENSIONS ntax=3;
       TAXLABELS
           A
           B
           C
       ;
   END;

   BEGIN Distances;
       DIMENSIONS ntax=3;
       FORMAT triangle=LOWER diagonal LABELS;
       MATRIX
       A 0.000000
       B 1.000000 0.000000
       C 2.000000 1.000000 0.000000
       ;
   END;

Example usage:

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np
   
   # Create distance matrix
   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   
   # Save to NEXUS (default format)
   dm.save("distances.nexus")
   
   # Load from NEXUS
   dm2 = DistanceMatrix.load("distances.nexus")
   
   # Export with upper triangular format
   nexus_str = dm.to_string(format="nexus", triangle="UPPER")

**MSA NEXUS Format**

The NEXUS format for sequence alignments includes a TAXA block and a CHARACTERS block:

.. code-block:: text

   #NEXUS

   BEGIN TAXA;
       DIMENSIONS ntax=2;
       TAXLABELS
           taxon1
           taxon2
       ;
   END;

   BEGIN CHARACTERS;
       DIMENSIONS nchar=8;
       FORMAT datatype=DNA missing=N gap=-;
       MATRIX
       taxon1    ACGTACGT
       taxon2    TGCAACGT
       ;
   END;

Example usage:

.. code-block:: python

   from phylozoo import MSA
   
   # Create MSA
   sequences = {"taxon1": "ACGTACGT", "taxon2": "TGCAACGT"}
   msa = MSA(sequences)
   
   # Save to NEXUS
   msa.save("alignment.nexus", format="nexus")
   
   # Load from NEXUS
   msa2 = MSA.load("alignment.nexus")

**Split System NEXUS Format**

The NEXUS format for split systems includes a TAXA block and a SPLITS block:

.. code-block:: text

   #NEXUS

   BEGIN TAXA;
       DIMENSIONS ntax=4;
       TAXLABELS
           A
           B
           C
           D
       ;
   END;

   BEGIN SPLITS;
       DIMENSIONS ntax=4;
       FORMAT labels=yes weights=yes;
       MATRIX
           1.0  A B | C D
           0.8  A C | B D
       ;
   END;

Example usage:

.. code-block:: python

   from phylozoo.core.split import SplitSystem, Split
   
   # Create split system
   split1 = Split({'A', 'B'}, {'C', 'D'})
   split2 = Split({'A', 'C'}, {'B', 'D'})
   splits = SplitSystem([split1, split2])
   
   # Save to NEXUS (default format)
   splits.save("splits.nexus")
   
   # Load from NEXUS
   splits2 = SplitSystem.load("splits.nexus")

FASTA Format
------------

FASTA is a simple text-based format for sequence data, widely used in bioinformatics.

**Classes using FASTA format:**

* **MSA** (default format)

**File Extensions:** ``.fasta``, ``.fa``, ``.fas``

**Format Structure**

Each sequence consists of a header line (starting with ``>``) followed by one or more 
lines of sequence data:

.. code-block:: text

   >taxon1
   ACGTACGT
   >taxon2
   TGCAACGT
   >taxon3
   AAAAACGT

Example usage:

.. code-block:: python

   from phylozoo import MSA
   
   # Create MSA
   sequences = {
       "taxon1": "ACGTACGT",
       "taxon2": "TGCAACGT",
       "taxon3": "AAAAACGT"
   }
   msa = MSA(sequences)
   
   # Save to FASTA (default format)
   msa.save("alignment.fasta")
   
   # Load from FASTA
   msa2 = MSA.load("alignment.fasta")
   
   # Customize line length (default is 80 characters)
   fasta_str = msa.to_string(format="fasta", line_length=60)

PHYLIP Format
-------------

PHYLIP is a compact format for distance matrices, compatible with standard PHYLIP tools.

**Classes using PHYLIP format:**

* **DistanceMatrix**

**File Extensions:** ``.phy``, ``.phylip``

**Format Structure**

The first line contains the number of taxa. Each subsequent line contains a taxon name 
(padded to 10 characters) followed by all distances for that taxon:

.. code-block:: text

   3
   A          0.00000 1.00000 2.00000
   B          1.00000 0.00000 1.00000
   C          2.00000 1.00000 0.00000

Example usage:

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np
   
   # Create distance matrix
   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   
   # Save to PHYLIP
   dm.save("distances.phy", format="phylip")
   
   # Load from PHYLIP
   dm2 = DistanceMatrix.load("distances.phy", format="phylip")

CSV Format
----------

CSV (Comma-Separated Values) is a simple tabular format, easy to read and edit in 
spreadsheet applications.

**Classes using CSV format:**

* **DistanceMatrix**

**File Extensions:** ``.csv``

**Format Structure**

The first row is a header with an empty first cell, followed by taxon labels. Each 
subsequent row contains a taxon label in the first column, followed by distances:

.. code-block:: text

   ,A,B,C
   A,0.0,1.0,2.0
   B,1.0,0.0,1.0
   C,2.0,1.0,0.0

Example usage:

.. code-block:: python

   from phylozoo import DistanceMatrix
   import numpy as np
   
   # Create distance matrix
   matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
   dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
   
   # Save to CSV
   dm.save("distances.csv", format="csv")
   
   # Load from CSV
   dm2 = DistanceMatrix.load("distances.csv", format="csv")
   
   # Use tab-separated values
   dm.save("distances.tsv", format="csv", delimiter="\t")

eNewick Format
--------------

eNewick (Extended Newick) is an extension of the standard Newick format that supports 
hybrid nodes and network structures with additional attributes.

**Classes using eNewick format:**

* **DirectedPhyNetwork** (default format)

**File Extensions:** ``.enewick``, ``.eNewick``, ``.enwk``, ``.nwk``, ``.newick``

**Format Structure**

eNewick extends Newick format to support hybrid nodes and edge attributes. It uses 
square brackets for hybrid node annotations and supports branch lengths, bootstrap 
values, and gamma probabilities:

.. code-block:: text

   ((A:1.0[&gamma=0.6],B:1.0[&gamma=0.4])#H1:0.5,C:2.0);

This example shows a hybrid node ``#H1`` with two incoming edges having gamma values 
0.6 and 0.4, and a branch length of 0.5.

Example usage:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Create network with hybrid node
   network = DirectedPhyNetwork(
       edges=[
           ("root", "u1"), ("root", "u2"),
           ("u1", "h", {"gamma": 0.6}),
           ("u2", "h", {"gamma": 0.4}),
           ("h", "leaf1")
       ],
       nodes=[("leaf1", {"label": "A"})]
   )
   
   # Save to eNewick (default format)
   network.save("network.enewick")
   
   # Load from eNewick
   network2 = DirectedPhyNetwork.load("network.enewick")

Newick Format
-------------

Newick is the standard format for phylogenetic trees and networks. PhyloZoo uses it 
for semi-directed networks.

**Classes using Newick format:**

* **SemiDirectedPhyNetwork** (default format)

**File Extensions:** ``.nwk``, ``.newick``, ``.enewick``, ``.eNewick``, ``.enw``

**Format Structure**

Standard Newick format for trees and networks:

.. code-block:: text

   ((A,B),C);

For semi-directed networks, undirected edges are represented in a way compatible with 
standard Newick parsers.

Example usage:

.. code-block:: python

   from phylozoo import SemiDirectedPhyNetwork
   
   # Create semi-directed network
   network = SemiDirectedPhyNetwork(
       undirected_edges=[(3, 1), (3, 2), (3, 4)],
       nodes=[
           (1, {"label": "A"}),
           (2, {"label": "B"}),
           (4, {"label": "C"})
       ]
   )
   
   # Save to Newick (default format)
   network.save("network.newick")
   
   # Load from Newick
   network2 = SemiDirectedPhyNetwork.load("network.newick")

DOT Format
----------

DOT is the Graphviz format for graph visualization. PhyloZoo supports DOT for network 
structures.

**Classes using DOT format:**

* **DirectedPhyNetwork**
* **DirectedMultiGraph** (default format)

**File Extensions:** ``.dot``, ``.gv``

**Format Structure**

DOT format uses a declarative syntax for graphs:

.. code-block:: text

   digraph {
       root -> u1;
       root -> u2;
       u1 -> h [label="gamma=0.6"];
       u2 -> h [label="gamma=0.4"];
       h -> leaf1 [label="A"];
   }

Example usage:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Create network
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B")],
       nodes=[("A", {"label": "A"}), ("B", {"label": "B"})]
   )
   
   # Save to DOT
   network.save("network.dot", format="dot")
   
   # Load from DOT
   network2 = DirectedPhyNetwork.load("network.dot", format="dot")

PhyloZoo-DOT Format
-------------------

PhyloZoo-DOT is a custom DOT format specifically designed for semi-directed networks, 
preserving the distinction between directed and undirected edges.

**Classes using PhyloZoo-DOT format:**

* **SemiDirectedPhyNetwork**
* **MixedMultiGraph** (default format)

**File Extensions:** ``.pzdot``

**Format Structure**

Similar to standard DOT but with special handling for undirected edges in semi-directed 
networks:

.. code-block:: text

   graph {
       node1 -- node2 [dir=none];
       node2 -> node3;
   }

Example usage:

.. code-block:: python

   from phylozoo import SemiDirectedPhyNetwork
   
   # Create semi-directed network
   network = SemiDirectedPhyNetwork(
       directed_edges=[(5, 4, {"gamma": 0.6})],
       undirected_edges=[(4, 1), (4, 2)],
       nodes=[
           (1, {"label": "A"}),
           (2, {"label": "B"})
       ]
   )
   
   # Save to PhyloZoo-DOT
   network.save("network.pzdot", format="phylozoo-dot")
   
   # Load from PhyloZoo-DOT
   network2 = SemiDirectedPhyNetwork.load("network.pzdot", format="phylozoo-dot")

PZ Format
---------

PZ is a custom JSON-based format for Squirrel quartet profile sets.

**Classes using PZ format:**

* **SqQuartetProfileSet** (default format)

**File Extensions:** ``.pz``

**Format Structure**

JSON-based format storing quartet profile sets with metadata:

.. code-block:: json

   {
       "profiles": [
           {
               "taxa": ["A", "B", "C", "D"],
               "quartets": [
                   {"split": ["A", "B"], "weight": 0.7},
                   {"split": ["A", "C"], "weight": 0.3}
               ]
           }
       ],
       "metadata": {...}
   }

Example usage:

.. code-block:: python

   from phylozoo.inference.squirrel import SqQuartetProfileSet, SqQuartetProfile
   from phylozoo.core.quartet import Quartet
   from phylozoo.core.split import Split
   
   # Create quartet profile set
   q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
   profile1 = SqQuartetProfile([q1])
   profileset = SqQuartetProfileSet(profiles=[profile1])
   
   # Save to PZ (default format)
   profileset.save("profiles.pz")
   
   # Load from PZ
   profileset2 = SqQuartetProfileSet.load("profiles.pz")

Edge List Format
----------------

Edge list is a simple text format for graphs, listing edges one per line.

**Classes using Edge List format:**

* **DirectedMultiGraph**

**File Extensions:** ``.edgelist``, ``.edges``

**Format Structure**

Each line contains an edge as two node identifiers:

.. code-block:: text

   1 2
   2 3
   3 1

Example usage:

.. code-block:: python

   from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
   
   # Create graph
   graph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (3, 1)])
   
   # Save to edge list
   graph.save("graph.edgelist", format="edgelist")
   
   # Load from edge list
   graph2 = DirectedMultiGraph.load("graph.edgelist", format="edgelist")

Error Handling
--------------

I/O operations raise appropriate exceptions when files cannot be read or written:

* ``PhyloZooIOError``: General I/O errors
* ``PhyloZooParseError``: Parsing errors (malformed files)
* ``PhyloZooFormatError``: Format-related errors (unsupported format, missing format handler)
* ``FileNotFoundError``: File not found errors

Example:

.. code-block:: python

   from phylozoo.utils.exceptions import PhyloZooParseError
   
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
   * :doc:`Networks (Basic) <../core/networks/directed/basic>` for network I/O
   * :doc:`Sequences <../core/sequences>` for MSA I/O
   * :doc:`Distance Matrices <../core/distance>` for distance matrix I/O
