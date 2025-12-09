# PhyloZoo Package Structure

This document describes the current package structure and organization.

## Overview

The package is organized into logical submodules based on functionality:

```
phylozoo/
├── __init__.py              # Main package initialization
├── structures/              # Fundamental data structures
│   ├── __init__.py
│   ├── splits.py           # Split, SplitSet, SplitSystem, QuartetSplitSet
│   ├── partition.py         # Partition class
│   └── circular.py         # CircularOrdering, CircularSetOrdering
│
├── diversity/               # Phylogenetic diversity calculations
│   ├── __init__.py
│   └── diversity.py        # DiversityCalculator, diversity functions
│
├── plotting/                # Network and tree visualization
│   ├── __init__.py
│   └── network_plot.py     # plot_network(), plot_tree()
│
├── squirrel/                # Network inference algorithms
│   ├── __init__.py
│   └── inference.py        # NetworkInferrer, inference functions
│
├── networks/               # Network types
│   ├── __init__.py
│   ├── directed.py         # DirectedNetwork
│   ├── semi_directed.py    # SemiDirectedNetwork, random_semi_directed_network
│   └── mixed_graph.py      # MixedGraph, MultiMixedGraph
│
├── utils/                   # Utility modules
│   ├── __init__.py
│   ├── _config.py          # Configuration (validation settings)
│   ├── tools.py            # General utility functions
│   └── distances.py        # DistanceMatrix
├── quarnet.py              # Quarnet classes (4-leaf networks)
├── quarnetset.py           # QuarnetSet classes
├── trinet.py               # Trinet classes (3-leaf networks)
├── trinetset.py            # TrinetSet classes
├── msa.py                  # Multiple Sequence Alignment
├── load_data.py            # Data loading functions
└── invariants/             # Data files directory
```

## Module Descriptions

### `structures/` - Fundamental Data Structures

**Purpose:** Contains fundamental data structures used throughout the package.

**Modules:**
- `splits.py`: Split, SplitSet, SplitSystem, QuartetSplitSet
- `partition.py`: Partition class
- `circular.py`: CircularOrdering, CircularSetOrdering

**Usage:**
```python
from phylozoo.structures import Split, SplitSystem, Partition, CircularOrdering
```

**Rationale:** These are fundamental concepts that other modules build upon. Grouping them makes dependencies clear and provides a clean foundation.

### `diversity/` - Phylogenetic Diversity

**Purpose:** Calculate phylogenetic diversity metrics.

**Modules:**
- `diversity.py`: DiversityCalculator class and diversity functions

**Usage:**
```python
from phylozoo.diversity import DiversityCalculator, phylogenetic_diversity
```

**Status:** Skeleton module - ready for implementation.

### `plotting/` - Visualization

**Purpose:** Plot and visualize phylogenetic networks and trees.

**Modules:**
- `network_plot.py`: plot_network(), plot_tree()
- `graph_plot.py`: plot_directed_multigraph(), plot_mixed_multigraph()

**Usage:**
```python
from phylozoo.visualize import plot_network, plot_tree
from phylozoo.visualize import plot_directed_multigraph, plot_mixed_multigraph
```

**Status:** Partially implemented - graph plotting functions available, network plotting placeholders.

### `squirrel/` - Network Inference

**Purpose:** Infer phylogenetic networks from data.

**Modules:**
- `inference.py`: NetworkInferrer class and inference functions

**Usage:**
```python
from phylozoo.squirrel import NetworkInferrer, infer_network_from_msa
```

**Status:** Skeleton module - ready for implementation.

### `networks/` - Network Types

**Purpose:** Classes for working with phylogenetic networks.

**Modules:**
- `directed.py`: DirectedNetwork
- `semi_directed.py`: SemiDirectedNetwork, random_semi_directed_network()
- `mixed_graph.py`: MixedGraph, MultiMixedGraph

**Usage:**
```python
from phylozoo.networks import DirectedNetwork, SemiDirectedNetwork, MixedGraph
```

### `utils/` - Utilities

**Purpose:** General utility functions and classes.

**Modules:**
- `_config.py`: Configuration settings (validation)
- `tools.py`: General utility functions
- `distances.py`: DistanceMatrix

## Import Strategy

### Recommended Imports

**Structures functionality:**
```python
from phylozoo.structures import Split, SplitSystem, Partition
```

**Networks:**
```python
from phylozoo.networks import DirectedNetwork, SemiDirectedNetwork
```

**Diversity:**
```python
from phylozoo.diversity import DiversityCalculator
```

**Plotting:**
```python
from phylozoo.visualize import plot_network
```

**Inference:**
```python
from phylozoo.squirrel import NetworkInferrer
```

### Backward Compatibility

To maintain backward compatibility with existing code, you may want to add re-exports in the main `__init__.py`:

```python
# In phylozoo/__init__.py
from .structures import Split, SplitSystem  # Re-export for convenience
```

## Future Organization Suggestions

Consider further organizing:

1. **Networks submodule:**
   - Move `dnetwork.py`, `sdnetwork.py` to `networks/`
   - Move `mixed_graph.py` from `utils/` to `networks/`

2. **Patterns submodule:**
   - Group `quarnet.py`, `quarnetset.py`, `trinet.py`, `trinetset.py` into `patterns/`

3. **IO submodule:**
   - Group `msa.py`, `load_data.py`, and invariant evaluation into `io/`

## Notes

- The `structures/` module contains fundamental structures that are used by all other modules
- New modules (`diversity/`, `plotting/`, `squirrel/`) are created as skeletons ready for implementation
- The structure is designed to be extensible and maintainable

