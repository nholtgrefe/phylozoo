Introduction
============

What is PhyloZoo?
-----------------

PhyloZoo is a Python package for working with phylogenetic networks. Phylogenetic networks 
extend phylogenetic trees by allowing for both divergence (splitting) and merging events, 
making them suitable for modeling processes like hybridization, horizontal gene transfer, 
and admixture :cite:`Huson2010`. Unlike traditional phylogenetic trees, networks can 
explicitly represent reticulate evolutionary events where lineages merge, providing a more 
realistic model of evolutionary history for many groups of organisms.

PhyloZoo provides a comprehensive toolkit for phylogenetic network analysis, including 
network representation, manipulation, inference from sequence data, diversity calculations, 
and visualization. The package is designed with a focus on correctness, performance, and 
ease of use.

Design Philosophy
-----------------

PhyloZoo follows several key design principles:

**Object-Oriented and Immutable**
   PhyloZoo uses an object-oriented design where core data structures (networks, quartets, 
   splits, etc.) are implemented as immutable classes. Once created, these objects cannot be 
   modified in-place, ensuring data integrity and making code more predictable and easier to 
   reason about. To modify a network, you create a new instance with the desired changes.

**Type Safety**
   All functions and classes use comprehensive type hints, enabling better IDE support, 
   static type checking with tools like mypy, and clearer documentation of expected inputs 
   and outputs.

**Comprehensive Documentation**
   All public functions, classes, and methods include NumPy-style docstrings with detailed 
   parameter descriptions, return values, exceptions, and examples. This ensures that the 
   code is self-documenting and accessible to users.

**Validation and Error Handling**
   PhyloZoo includes a custom exception hierarchy for clear, specific error messages. Network 
   validation ensures that objects always represent valid phylogenetic structures, catching 
   errors early and providing helpful diagnostic information.

**Performance**
   Where appropriate, PhyloZoo leverages NumPy for efficient numerical operations and 
   supports optional Numba JIT compilation for computationally intensive algorithms. The 
   package is designed to handle both small-scale exploratory analyses and larger 
   production workflows.

Package Structure
-----------------

PhyloZoo is organized into several main modules:

**Core Module** (`phylozoo.core`)
   The core module contains fundamental data structures and classes:
   
   * **Networks**: ``DirectedPhyNetwork`` and ``SemiDirectedPhyNetwork`` classes for 
     representing phylogenetic networks. Directed networks are fully directed DAGs with 
     explicit root and hybrid nodes, while semi-directed networks allow undirected tree 
     edges for more flexible representation :cite:`NetworkRepresentation2024`.
   
   * **Quartets**: ``Quartet``, ``QuartetProfile``, and ``QuartetProfileSet`` classes for 
     working with four-taxon relationships, which are fundamental building blocks for 
     network inference :cite:`QuartetInference2024`.
   
   * **Splits**: ``Split`` and ``SplitSystem`` classes for representing bipartitions of 
     taxa, a common way to encode phylogenetic relationships.
   
   * **Sequences**: ``MSA`` (Multiple Sequence Alignment) class with efficient NumPy-based 
     storage and bootstrapping capabilities.
   
   * **Distance**: ``DistanceMatrix`` class for pairwise distance data with support for 
     various distance matrix properties and classifications.
   
   * **Primitives**: Fundamental structures like ``Partition``, ``CircularOrdering``, and 
     ``CircularSetOrdering`` used throughout the package.

**Inference Module** (`phylozoo.inference`)
   Algorithms for inferring phylogenetic networks from data:
   
   * **Squirrel**: Implementation of the SQuaRE (SQuirrel quartet-based REconstruction) 
     algorithm for inferring semi-directed phylogenetic networks from quartet profiles 
     :cite:`SquirrelAlgorithm2024`. This includes quartet profile generation, cycle 
     resolution, and network assembly.
   
   * **Delta Heuristic**: Methods for generating quartet profiles from distance matrices 
     using the delta heuristic approach.

**Panda Module** (`phylozoo.panda`)
   Phylogenetic diversity calculations and optimization:
   
   * **Diversity Measures**: Protocol-based architecture for defining and computing 
     various phylogenetic diversity metrics, such as all-paths diversity.
   
   * **Optimization**: Algorithms for finding sets of taxa that maximize phylogenetic 
     diversity, including greedy approaches and exact solutions using dynamic programming 
     for certain network structures :cite:`DiversityOptimization2024`.

**Visualization Module** (`phylozoo.viz`)
   Flexible plotting system for networks:
   
   * **Network Plotting**: Functions for visualizing directed and semi-directed networks 
     with customizable layouts, styling, and multiple backend support (Matplotlib, PyQtGraph).
   
   * **Layout Algorithms**: Various layout algorithms including combining layouts for 
     networks with tree-of-blobs structure.

**Utilities Module** (`phylozoo.utils`)
   Supporting functionality:
   
   * **Exceptions**: Comprehensive custom exception hierarchy for clear error reporting.
   
   * **Validation**: Network validation utilities and decorators.
   
   * **I/O**: File format support including eNewick, DOT, and other formats.

Key Features
------------

* **Network Representation**: Support for both directed and semi-directed phylogenetic 
  networks with comprehensive validation and feature extraction.

* **Network Manipulation**: Tools for creating, modifying, and transforming networks, 
  including binary resolution, subnetwork extraction, and network generators.

* **Sequence Analysis**: Efficient MSA handling with NumPy-based storage, bootstrapping 
  support, and distance calculations.

* **Network Inference**: Implementation of quartet-based network inference algorithms, 
  including the SQuaRE method for reconstructing networks from quartet profiles.

* **Visualization**: Flexible plotting system with multiple layout algorithms and backend 
  support for creating publication-quality figures.

* **Diversity Analysis**: Phylogenetic diversity calculations with optimization algorithms 
  for selecting diverse sets of taxa.

* **Type Safety**: Comprehensive type hints throughout the codebase for better IDE support 
  and static type checking.

Quick Start
-----------

Here's a simple example to get started:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # Create a simple network
   network = DirectedPhyNetwork(
       edges=[("root", "A"), ("root", "B")],
       nodes=[
           ("A", {"label": "A"}),
           ("B", {"label": "B"})
       ]
   )
   
   print(f"Network has {network.num_nodes} nodes and {network.num_edges} edges")
   print(f"Leaves: {network.leaves}")

Core Concepts
-------------

**Directed Networks**: Fully directed phylogenetic networks with a single root, where 
hybrid nodes have in-degree >= 2. These networks explicitly represent the direction of 
evolutionary time and are suitable for rooted phylogenetic analyses.

**Semi-Directed Networks**: Networks with directed hybrid edges and undirected tree edges, 
allowing for more flexible representation of evolutionary relationships. These are 
particularly useful for unrooted analyses and can represent the same underlying 
phylogenetic structure in multiple ways.

**Quartets**: Four-taxon unrooted trees, fundamental building blocks for network inference. 
Quartet-based methods are widely used because they provide a natural way to decompose 
complex phylogenetic relationships into smaller, more manageable pieces.

**Splits**: Bipartitions of taxa, used to represent evolutionary relationships. A split 
divides the set of taxa into two non-empty groups, representing a potential evolutionary 
divergence event.

**Blobs**: Maximal biconnected components of a network. Blobs represent regions of 
reticulation and are important for understanding network structure and complexity.

For more detailed information, see the :doc:`Library <../library/index>` section.
