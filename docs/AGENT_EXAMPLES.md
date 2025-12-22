# AGENT_EXAMPLES.md — Concrete examples and templates

This file contains copy-pasteable examples that agents should use as templates.

1) DistanceMatrix usage (immutable container)

```python
from phylozoo.core.distance.base import DistanceMatrix
import numpy as np

matrix = np.array([[0, 1], [1, 0]], dtype=float)
dm = DistanceMatrix(matrix, labels=['A', 'B'])
assert dm.np_array.flags.writeable is False
assert dm.get_distance('A', 'B') == 1.0
```

Note: `DistanceMatrix` builds an internal dict for O(1) lookups. Labels must be hashable (e.g. int, str, tuple). If you pass mutable/unhashable labels (list, set, dict) canonicalize them first or convert to a hashable type.

2) NumPy/Numba safe function pattern

- Keep Numba-compiled functions simple, with fixed dtypes and contiguous arrays.
- Write a pure-Python fallback for readability and tests.

Example pattern:

```python
import numpy as np
from numba import njit

@njit
def _fast_dot(a: np.ndarray, b: np.ndarray) -> float:
    s = 0.0
    for i in range(a.shape[0]):
        s += a[i] * b[i]
    return s

# Pure-Python wrapper used in tests and as fallback
def fast_dot(a: np.ndarray, b: np.ndarray) -> float:
    a = np.ascontiguousarray(a, dtype=np.float64)
    b = np.ascontiguousarray(b, dtype=np.float64)
    return float(_fast_dot(a, b))
```

3) Docstring template (NumPy style)

Use this in new functions and classes (examples in tests should mirror these):

```python
def add(a: int, b: int) -> int:
    """
    Add two integers.

    Parameters
    ----------
    a : int
        First integer.
    b : int
        Second integer.

    Returns
    -------
    int
        The sum of `a` and `b`.

    Examples
    --------
    >>> add(1, 2)
    3
    """
    return a + b
```

4) Pytest test template

```python
import pytest
from phylozoo.core.distance.base import DistanceMatrix
import numpy as np

def test_distance_matrix_basic() -> None:
    """Basic behavior: construction and symmetry."""
    m = np.array([[0, 2], [2, 0]], dtype=float)
    dm = DistanceMatrix(m, labels=['X','Y'])
    assert len(dm) == 2
    assert dm.get_distance('X', 'Y') == 2.0
```

5) Running individual tests and markers

- Run a single test module:

```bash
pytest tests/test_splits.py -q
```

- Run tests excluding slow ones:

```bash
pytest -q -m "not slow"
```

6) Example: micro-benchmark pattern

- Provide a small timed comparison in `sandbox/findings.txt` or as a test marked `slow` to demonstrate performance parity.

```python
import time
start = time.perf_counter()
# run function
elapsed = time.perf_counter() - start
print(f"Elapsed: {elapsed:.6f}s")
```

Use these templates when adding features or fixing bugs; they reduce review friction.