# Examples

This directory contains example scripts demonstrating how to use the phylozoo package.

## Running Examples

To run an example:

```bash
python examples/example_basic.py
```

Or from the project root:

```bash
python -m examples.example_basic
```

## Available Examples

- `example_basic.py` - Basic usage examples showing how to create and work with splits, networks, and MSAs
- `example_directed_phynetwork.py` - Comprehensive examples demonstrating DirectedPhyNetwork features including trees, hybrid networks, edge attributes (branch_length, bootstrap, gamma), and network topology properties

## Creating New Examples

When creating new examples:

1. Use descriptive filenames (e.g., `example_<feature>.py`)
2. Include a docstring at the top explaining what the example demonstrates
3. Use type hints for all functions
4. Include a `main()` function that can be called directly
5. Add a `if __name__ == "__main__":` block to allow direct execution

Example structure:

```python
"""
Example demonstrating [feature].

This example shows how to [do something].
"""
from phylozoo import SomeClass


def main() -> None:
    """Run the example."""
    # Example code here
    pass


if __name__ == "__main__":
    main()
```

