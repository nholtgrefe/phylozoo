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

Why PhyloZoo?
-------------

PhyloZoo offers several advantages for phylogenetic network analysis:

**Unified Framework**
   PhyloZoo provides a single, consistent interface for working with networks, quartets, 
   splits, sequences, and distance matrices. All data structures follow the same design 
   principles, making it easy to move between different representations and combine them 
   in your analyses.

**Robust Validation**
   All network objects are automatically validated upon creation, ensuring that you're 
   always working with valid phylogenetic structures. This catches errors early and 
   provides clear, helpful error messages when something goes wrong.

**Comprehensive Functionality**
   From network representation and manipulation to visualization, PhyloZoo provides 
   support for phylogenetic network workflows. You can work with sequence data, 
   networks, quartets, and splits without switching between different tools.

**Performance and Scalability**
   PhyloZoo is designed to handle both small exploratory analyses and larger production 
   workflows. Efficient NumPy-based operations and optional parallelization support make 
   it suitable for analyzing networks of various sizes.

**Extensible Design**
   The package uses a protocol-based architecture for diversity measures and provides 
   flexible layout and styling systems for visualization, making it easy to extend and 
   customize for your specific needs.

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
     edges for more flexible representation :cite:`NetworkRepresentation2024`. 
     See :doc:`Networks <core/networks/index>` for details.
   
   * **Quartets**: ``Quartet``, ``QuartetProfile``, and ``QuartetProfileSet`` classes for 
     working with four-taxon relationships, which are fundamental building blocks for 
     network inference :cite:`QuartetInference2024`. See :doc:`Quartets <core/quartets/index>` 
     for details.
   
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

**Visualization Module** (`phylozoo.viz`)
   Flexible plotting system for networks. See :doc:`Visualization Module <visualization/index>` 
   for detailed documentation.
   
   * **Network Plotting**: Functions for visualizing directed and semi-directed networks 
     with customizable layouts and styling using Matplotlib. See :doc:`Plotting <visualization/plotting>` 
     for details.
   
   * **Layout Algorithms**: Custom PhyloZoo layouts (pz-dag, pz-radial) and access to 
     standard NetworkX and Graphviz layouts for various visualization needs. See 
     :doc:`Layouts <visualization/layouts>` for details.

**Utilities Module** (`phylozoo.utils`)
   Supporting functionality. See :doc:`Utils Module <utils/index>` for detailed documentation.
   
   * **Exceptions**: Comprehensive custom exception hierarchy for clear error reporting. 
     See :doc:`Exceptions <utils/exceptions>` for details.
   
   * **Validation**: Network validation utilities and decorators. See 
     :doc:`Validation <utils/validation>` for details.
   
   * **I/O**: File format support including eNewick, DOT, and other formats. See 
     :doc:`I/O <utils/io>` for details.
   
   * **Parallelization**: Support for parallel execution of computationally intensive tasks. 
     See :doc:`Parallelization <utils/parallelization>` for details.

Key Features
------------

* **Network Representation**: Support for both directed and semi-directed phylogenetic 
  networks with comprehensive validation and feature extraction. Networks can be created 
  programmatically or loaded from various file formats.

* **Network Manipulation**: Tools for creating, modifying, and transforming networks, 
  including binary resolution, subnetwork extraction, and network generators for creating 
  test networks.

* **Sequence Analysis**: Efficient MSA (Multiple Sequence Alignment) handling with 
  NumPy-based storage, bootstrapping support, and distance calculations. Supports 
  common formats like FASTA and NEXUS.

* **Visualization**: Flexible plotting system with multiple layout algorithms (including 
  custom PhyloZoo layouts and standard NetworkX/Graphviz layouts) for creating 
  publication-quality figures. Supports customizable styling and rendering.

* **Parallelization**: Built-in support for parallel execution of computationally 
  intensive tasks, allowing you to leverage multiple CPU cores for faster analysis.

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

To get started with PhyloZoo, see the :doc:`Quickstart <../quickstart>` guide for installation 
instructions and basic examples. For detailed documentation on specific modules, see:

* :doc:`Core Module <core/index>`: Data structures and network operations
* :doc:`Visualization Module <visualization/index>`: Network plotting and visualization

For complete API reference, see the :doc:`API Reference <../api/index>` section.
