Package Structure
================

This document describes the package structure, import strategy, and design principles of PhyloZoo.

Package Organization
--------------------

PhyloZoo is organized into logical modules based on functionality:

**Core Module** (`phylozoo.core`)
   Fundamental data structures and classes:
   
   * **Networks** (`phylozoo.core.network`): ``DirectedPhyNetwork`` and ``SemiDirectedPhyNetwork``
   * **Quartets** (`phylozoo.core.quartet`): ``Quartet``, ``QuartetProfile``, ``QuartetProfileSet``
   * **Splits** (`phylozoo.core.split`): ``Split``, ``SplitSystem``
   * **Sequences** (`phylozoo.core.sequence`): ``MSA`` (Multiple Sequence Alignment)
   * **Distance** (`phylozoo.core.distance`): ``DistanceMatrix``
   * **Primitives** (`phylozoo.core.primitives`): ``Partition``, ``CircularOrdering``, ``CircularSetOrdering``

**Inference Module** (`phylozoo.inference`)
   Network inference algorithms:
   
   * **Squirrel** (`phylozoo.inference.squirrel`): Squirrel algorithm for network reconstruction

**Diversity Module** (`phylozoo.panda`)
   Phylogenetic diversity calculations:
   
   * **Diversity Measures**: Protocol-based architecture for diversity metrics

**Visualization Module** (`phylozoo.viz`)
   Network plotting and visualization:
   
   * **Network Plotting**: Functions for directed and semi-directed networks
   * **Layout Algorithms**: Custom and standard layouts
   * **Styling**: Customizable visual appearance

**Utilities Module** (`phylozoo.utils`)
   Supporting functionality:
   
   * **Exceptions**: Custom exception hierarchy
   * **Validation**: Network validation utilities
   * **I/O**: File format support
   * **Parallelization**: Parallel execution support

Import Strategy
---------------

PhyloZoo follows a consistent import strategy:

**Main Datatypes (Top-Level)**
   Core data structures are available at the package top-level:
   
   .. code-block:: python
   
      from phylozoo import (
          DirectedPhyNetwork,
          SemiDirectedPhyNetwork,
          Quartet,
          QuartetProfile,
          QuartetProfileSet,
          Split,
          SplitSystem,
          DistanceMatrix,
          MSA,
      )

**Module-Level Imports**
   All public functions and classes are available at their folder level:
   
   .. code-block:: python
   
      # Core modules
      from phylozoo.core.quartet import Quartet, QuartetProfile
      from phylozoo.core.sequence import bootstrap, hamming_distances
      from phylozoo.core.distance import is_kalmanson
      
      # Inference
      from phylozoo.inference.squirrel import squirrel, delta_heuristic
      
      # Diversity
      from phylozoo.panda import diversity, AllPathsDiversity
      
      # Visualization
      from phylozoo.viz import plot_dnetwork, plot_sdnetwork

**Network Submodules (Aliases)**
   Network submodules have convenient aliases to avoid deep nesting:
   
   .. code-block:: python
   
      # Instead of: from phylozoo.core.network.dnetwork import ...
      from phylozoo.core.dnetwork import DirectedPhyNetwork, classifications, features
      
      # Instead of: from phylozoo.core.network.sdnetwork import ...
      from phylozoo.core.sdnetwork import SemiDirectedPhyNetwork, classifications, features

**Exceptions (Top-Level)**
   All custom exceptions are available at the package top-level:
   
   .. code-block:: python
   
      from phylozoo import (
          PhyloZooError,
          PhyloZooValueError,
          PhyloZooNetworkError,
          # ... other exceptions
      )

**Internal Items**
   Items starting with `_` are considered internal and are **not** re-exported at higher levels. 
   They can only be imported directly from their source file.

Design Principles
-----------------

PhyloZoo follows several key design principles:

**Immutability**
   Core data structures (networks, quartets, splits, etc.) are immutable. Once created, 
   objects cannot be modified in-place. To modify a network, you create a new instance 
   with the desired changes. This ensures data integrity and makes code more predictable.

**Type Safety**
   All functions and classes use comprehensive type hints, enabling better IDE support, 
   static type checking with tools like mypy, and clearer documentation of expected inputs 
   and outputs.

**Validation**
   All network objects are automatically validated upon creation, ensuring that you're 
   always working with valid phylogenetic structures. This catches errors early and 
   provides clear, helpful error messages.

**Comprehensive Documentation**
   All public functions, classes, and methods include NumPy-style docstrings with detailed 
   parameter descriptions, return values, exceptions, and examples.

**Performance**
   PhyloZoo leverages NumPy for efficient numerical operations and supports optional 
   Numba JIT compilation for computationally intensive algorithms. The package is designed 
   to handle both small-scale exploratory analyses and larger production workflows.

**Modular Architecture**
   The package uses a modular architecture with clear separation of concerns. Layout 
   computation, styling, and rendering are separated in the visualization module, and 
   diversity measures use a protocol-based architecture for extensibility.

**No Deep Imports**
   Do not import from deeper levels than submodules (e.g., don't import from `backends` or 
   other internal subdirectories). Use the submodule level instead.

**No Package-Level Namespace Shortcuts**
   We intentionally do **not** provide convenience module namespaces like `phylozoo.sequence` 
   or `from phylozoo import sequence`. If you want a module namespace, import it directly 
   from its location, e.g., `from phylozoo.core import sequence`.

For more details on the import strategy, see the :doc:`Import Strategy <import_strategy>` 
documentation (if available).
