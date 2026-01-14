# PhyloZoo Import Strategy

This document describes the import strategy and API structure for the PhyloZoo package.

## General Principles

1. **All public functions and classes are available at their folder level**
   - Example: `from phylozoo.inference.squirrel import SqQuartetProfile`
   - Internal items (names starting with `_`) are **not** re-exported to higher levels

2. **Main datatypes are available at package top-level**
   - Core data structures can be imported directly: `from phylozoo import Quartet, Split, DistanceMatrix, MSA`

3. **Exceptions are available at package top-level**
   - All custom exceptions: `from phylozoo import PhyloZooError, PhyloZooValueError, ...`

4. **No package-level namespace shortcuts**
   - We intentionally do **not** provide convenience module namespaces like `phylozoo.sequence` or `from phylozoo import sequence`.
   - If you want a module namespace, import it directly from its location, e.g. `from phylozoo.core import sequence` or `from phylozoo.inference import squirrel`.

## Module-Specific Import Patterns

### Core Module (`phylozoo.core`)

#### Main Datatypes (Top-Level)
The main data structures are re-exported at the package top-level:

```python
from phylozoo import (
    # Networks
    DirectedPhyNetwork,
    SemiDirectedPhyNetwork,
    # Structures
    Split,
    SplitSystem,
    Partition,
    CircularOrdering,
    CircularSetOrdering,
    # Distance
    DistanceMatrix,
    # Sequence
    MSA,
    # Quartet
    Quartet,
    QuartetProfile,
    QuartetProfileSet,
)
```

#### Network Submodules (Aliases)
The network submodules have convenient aliases to avoid deep nesting:

```python
# Instead of: from phylozoo.core.network.dnetwork import ...
from phylozoo.core.dnetwork import DirectedPhyNetwork, classifications, features

# Instead of: from phylozoo.core.network.sdnetwork import ...
from phylozoo.core.sdnetwork import SemiDirectedPhyNetwork, classifications, features
```

#### Other Core Modules
All other core modules are available at their folder level:

```python
# Quartet module
from phylozoo.core.quartet import Quartet, QuartetProfile, QuartetProfileSet

# Sequence module
from phylozoo.core.sequence import MSA, bootstrap, bootstrap_per_gene, hamming_distances

# Distance module
from phylozoo.core.distance import DistanceMatrix, classifications, operations

# Split module
from phylozoo.core.split import Split, SplitSystem, WeightedSplitSystem

# Primitives module
from phylozoo.core.primitives import Partition, CircularOrdering, CircularSetOrdering
```

**Note:** We do not re-export these core subpackages at the package level (i.e., `phylozoo.sequence` is not part of the public API).

### Inference Module (`phylozoo.inference`)

Everything is available at the folder level:

```python
# Squirrel submodule
from phylozoo.inference.squirrel import (
    SqQuartetProfile,
    SqQuartetProfileSet,
    adapted_quartet_joining,
    quartet_joining,
    resolve_cycles,
    sqprofileset_similarity,
    sqprofileset_from_network,
    delta_heuristic,
    squirrel,
    # ... other public functions
)

# Utils submodule
from phylozoo.inference.utils import root_at_outgroup
```

**Note:** Internal functions (starting with `_`) are **not** exported and cannot be imported from the module level.

### Panda Module (`phylozoo.panda`)

Everything is available directly from `phylozoo.panda` (no need to go into `measure` subfolder):

```python
from phylozoo.panda import (
    diversity,
    greedy_max_diversity,
    marginal_diversities,
    solve_max_diversity,
    all_paths,
    AllPathsDiversity,
    DiversityMeasure,
)
```

### Utils Module (`phylozoo.utils`)

The utils module is **not** re-exported at the package level. Access utilities directly:

```python
from phylozoo.utils import exceptions, validation, io
from phylozoo.utils.exceptions import PhyloZooError
from phylozoo.utils.validation import validate_network
```

### Exceptions (`phylozoo.utils.exceptions`)

All exceptions are available at the package top-level:

```python
from phylozoo import (
    PhyloZooError,
    PhyloZooNotImplementedError,
    PhyloZooValueError,
    PhyloZooTypeError,
    PhyloZooRuntimeError,
    PhyloZooImportError,
    PhyloZooAttributeError,
    PhyloZooNetworkError,
    PhyloZooNetworkStructureError,
    PhyloZooNetworkDegreeError,
    PhyloZooNetworkAttributeError,
    PhyloZooIOError,
    PhyloZooParseError,
    PhyloZooFormatError,
    PhyloZooAlgorithmError,
    PhyloZooVisualizationError,
    PhyloZooLayoutError,
    PhyloZooBackendError,
    PhyloZooStateError,
    PhyloZooWarning,
    PhyloZooIdentifierWarning,
    PhyloZooEmptyNetworkWarning,
    PhyloZooSingleNodeNetworkWarning,
    PhyloZooGeneratorError,
    PhyloZooGeneratorStructureError,
    PhyloZooGeneratorDegreeError,
)
```

### Viz Module (`phylozoo.viz`)

Everything is available at the submodule level (not from deeper levels like `backends`):

```python
# Directed network plotting
from phylozoo.viz.dnetwork import (
    plot_network,
    plot_tree,
    compute_dag_layout,
    NetworkStyle,
    default_style,
    get_backend,
    register_backend,
)

# Semi-directed network plotting
from phylozoo.viz.sdnetwork import ...

# Graph plotting
from phylozoo.viz.graphs import plot_directed_multigraph, plot_mixed_multigraph
```

**Note:** Do not import from deeper levels like `phylozoo.viz.dnetwork.backends` - use the submodule level instead.

## Internal Items

Items starting with `_` are considered internal and are **not** re-exported at higher levels. They can only be imported directly from their source file:

```python
# This will NOT work (internal item):
from phylozoo.inference.squirrel import _qprofiles_to_circular_ordering  # ImportError

# This WILL work (if you really need it):
from phylozoo.inference.squirrel.cycle_resolution import _qprofiles_to_circular_ordering
```

## Examples

### Complete Example

```python
# Top-level imports (main datatypes)
from phylozoo import (
    DirectedPhyNetwork,
    SemiDirectedPhyNetwork,
    Quartet,
    Split,
    DistanceMatrix,
    MSA,
    PhyloZooError,
    PhyloZooValueError,
)

# Core module imports
from phylozoo.core.dnetwork import classifications, features
from phylozoo.core.sequence import bootstrap, bootstrap_per_gene
from phylozoo.core.distance import is_kalmanson

# Inference module
from phylozoo.inference.squirrel import SqQuartetProfile, squirrel

# Panda module
from phylozoo.panda import diversity, AllPathsDiversity

# Viz module
from phylozoo.viz.dnetwork import plot_network, NetworkStyle
```

## Summary

- **Top-level**: Main datatypes and exceptions (`from phylozoo import ...`)
- **Module level**: All public functions/classes (`from phylozoo.module import ...`)
- **Aliases**: `phylozoo.core.dnetwork` and `phylozoo.core.sdnetwork` for convenience
- **No deep imports**: Don't go deeper than submodule level (e.g., don't import from `backends`)
- **No internal exports**: Items starting with `_` are not re-exported
