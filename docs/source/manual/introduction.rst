Introduction
============

What is PhyloZoo?
-----------------

PhyloZoo is a Python package for working with phylogenetic networks. Phylogenetic networks 
extend phylogenetic trees by allowing for both divergence (splitting) and merging events, 
making them suitable for modeling processes like hybridization, horizontal gene transfer, 
and admixture. Unlike traditional phylogenetic trees, networks can 
explicitly represent reticulate evolutionary events where lineages merge, providing a more 
realistic model of evolutionary history for many groups of organisms.

PhyloZoo provides a comprehensive toolkit for phylogenetic network analysis, including 
network representation, manipulation, and visualization. The package is designed with a focus on correctness, performance, and 
ease of use. It supports both directed and semi-directed network representations, allowing 
you to work with rooted and unrooted phylogenetic analyses. PhyloZoo integrates seamlessly 
with the Python scientific computing ecosystem, using NumPy for efficient numerical 
operations and providing a clean, intuitive API for phylogenetic analysis workflows.

Beyond network representation, PhyloZoo offers native support for quartets, split systems,
multiple sequence alignments, and distance matrices within a consistent interface.
Conversions between these representations are supported, allowing analyses to move
flexibly between data types as required. All core data structures are validated upon
construction to ensure well-defined phylogenetic objects, improving reliability and
reproducibility. The package includes support for standard file formats such as eNewick,
DOT, FASTA, and NEXUS, and flexible visualization functionality with
customizable layouts for figures.

Package Structure
-----------------

PhyloZoo is organized into several main modules. For detailed documentation on each module, 
see the corresponding sections in this manual:

**Core Module** (`phylozoo.core`)
   The core module contains fundamental data structures and classes. See :doc:`Core Module <core/index>` 
   for detailed documentation.
   
   * **Networks**: ``DirectedPhyNetwork`` and ``SemiDirectedPhyNetwork`` classes for 
     representing phylogenetic networks. Directed networks are fully directed DAGs with 
     explicit root and hybrid nodes, while semi-directed networks allow undirected tree 
     edges for modelling root uncertainty :cite:`NetworkRepresentation2024`. 
     See :doc:`Networks <core/networks/index>` for details.
   
   * **Quartets**: ``Quartet``, ``QuartetProfile``, and ``QuartetProfileSet`` classes for 
     working with four-taxon relationships, which are fundamental building blocks for 
     network inference. See :doc:`Quartets <core/quartets/index>` for details.
   
   * **Splits**: ``Split`` and ``SplitSystem`` classes for representing bipartitions of 
     taxa, a common way to encode phylogenetic relationships. See :doc:`Splits <core/splits/index>` 
     for details.
   
   * **Sequences**: ``MSA`` (Multiple Sequence Alignment) class with efficient NumPy-based 
     storage and bootstrapping capabilities. See :doc:`Sequences <core/sequences>` for details.
   
   * **Distance**: ``DistanceMatrix`` class for pairwise distance data with support for 
     various distance matrix properties and classifications. See :doc:`Distance <core/distance>` 
     for details.
   
   * **Primitives**: Fundamental structures like ``Partition``, ``CircularOrdering``, and 
     ``CircularSetOrdering`` used throughout the package. See :doc:`Primitives <core/primitives/index>` 
     for details.

**Visualization Module** (:mod:`phylozoo.viz`)
   Flexible plotting system for networks. See :doc:`Visualization Module <visualization/index>` 
   for detailed documentation.
   
   * **Network Plotting**: Functions for visualizing directed and semi-directed networks 
     with customizable layouts and styling using Matplotlib. See :doc:`Plotting <visualization/plotting>` 
     for details.
   
   * **Layout Algorithms**: Custom PhyloZoo layouts (pz-dag, pz-radial) and access to 
     standard NetworkX and Graphviz layouts for various visualization needs. See 
     :doc:`Layouts <visualization/layouts>` for details.

**Utilities Module** (:mod:`phylozoo.utils`)
   Supporting functionality. See :doc:`Utils Module <utils/index>` for detailed documentation.
   
   * **Exceptions**: Comprehensive custom exception hierarchy for clear error reporting. 
     See :doc:`Exceptions <utils/exceptions>` for details.
   
   * **Validation**: Class and object validation utilities and decorators. See 
     :doc:`Validation <utils/validation>` for details.
   
   * **I/O**: File format support including eNewick, DOT, and many other formats. See 
     :doc:`I/O <utils/io>` for details.

Design Philosophy
-----------------

PhyloZoo follows several key design principles that benefit end users:

**Object-Oriented and Immutable**
   PhyloZoo uses an object-oriented design where core data structures (networks, quartets, 
   splits, etc.) are implemented as immutable classes. Once created, these objects cannot be 
   modified in-place, ensuring data integrity and making code more predictable and easier to 
   reason about. To modify a network, you create a new instance with the desired changes.

**Comprehensive Documentation**
   All public functions, classes, and methods include detailed docstrings with parameter 
   descriptions, return values, exceptions, and examples. This ensures that the code is 
   self-documenting and accessible to users.

**Validation and Error Handling**
   PhyloZoo includes a custom exception hierarchy for clear, specific error messages. Network 
   validation ensures that objects always represent valid phylogenetic structures, catching 
   errors early and providing helpful diagnostic information.

**Performance**
   Where appropriate, PhyloZoo leverages NumPy for efficient numerical operations and 
   supports optional Numba JIT compilation for computationally intensive algorithms. The 
   package is designed to handle both small-scale exploratory analyses and larger 
   production workflows.

Getting Started
---------------

To get started with PhyloZoo, see the :doc:`Installation Guide <installation>` for detailed installation 
instructions and visit the detailed documentation on specific modules:

* :doc:`Core Module <core/index>`: Data structures and network operations
* :doc:`Visualization Module <visualization/index>`: Network plotting and visualization
* :doc:`Utilities Module <utils/index>`: Utility functions and classes

Alternatively, see the :doc:`Quickstart <../tutorials/quickstart>` guide for a quickstart tutorial.

For complete API reference, see the :doc:`API Reference <../api/index>` section.
