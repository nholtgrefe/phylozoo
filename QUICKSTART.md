# PhyloZoo Quick Start Guide

This guide will help you get started with the phylozoo package.

## Initial Setup

### 1. Create a Virtual Environment (Recommended)

```bash
cd /home/nholtgreve/Documents/Code/phylozoo
python -m venv .venv
source .venv/bin/activate  # On Linux/Mac
# or
.venv\Scripts\activate  # On Windows
```

### 2. Install the Package

Install in development mode:

```bash
pip install -e .
```

For development with all tools (testing, linting, type checking):

```bash
pip install -e ".[dev]"
```

Or just for testing:

```bash
pip install -e ".[test]"
```

## Running Tests

### Run All Tests

```bash
pytest
```

### Run Tests with Coverage

```bash
pytest --cov=phylozoo --cov-report=html
```

This will generate an HTML coverage report in `htmlcov/index.html`.

### Run Specific Test Files

```bash
pytest tests/test_splits.py
pytest tests/test_utils.py
pytest tests/test_networks.py
```

### Run Tests Matching a Pattern

```bash
pytest -k "test_split"
```

### Run Tests in Verbose Mode

```bash
pytest -v
```

### Run Tests with Print Statements

```bash
pytest -s
```

## Running Examples

### Basic Example

```bash
python examples/example_basic.py
```

Or:

```bash
python -m examples.example_basic
```

## Running Scripts

```bash
python scripts/example_script.py -i input.txt -o output.txt
```

## Development Tools

### Type Checking with MyPy

```bash
mypy src/phylozoo
```

### Code Formatting with Black

```bash
black src/phylozoo tests examples scripts
```

### Linting with Ruff

```bash
ruff check src/phylozoo tests
```

## Project Structure

```
phylozoo/
├── src/
│   └── phylozoo/          # Main package
│       ├── __init__.py
│       ├── dnetwork.py
│       ├── invariant_eval.py
│       ├── load_data.py
│       ├── msa.py
│       ├── quarnet.py
│       ├── quarnetset.py
│       ├── sdnetwork.py
│       ├── splits.py
│       ├── trinet.py
│       ├── trinetset.py
│       ├── invariants/     # Data files directory
│       └── utils/          # Utility modules
│           ├── __init__.py
│           ├── _config.py
│           ├── circular.py
│           ├── distances.py
│           ├── mixed_graph.py
│           ├── partition.py
│           ├── polynomial.py
│           ├── polynomial2.py
│           └── tools.py
├── tests/                  # Test suite
│   ├── __init__.py
│   ├── conftest.py
│   ├── test_splits.py
│   ├── test_utils.py
│   └── test_networks.py
├── examples/               # Example scripts
│   ├── example_basic.py
│   └── README.md
├── scripts/                # Utility scripts
│   ├── example_script.py
│   └── README.md
├── pyproject.toml          # Package configuration
├── README.md               # Main documentation
├── QUICKSTART.md           # This file
└── .gitignore              # Git ignore rules
```

## Next Steps

1. **Implement your code**: Replace placeholder implementations in the modules with your actual code
2. **Add tests**: Write comprehensive tests for your implementations
3. **Add examples**: Create more example scripts demonstrating package features
4. **Documentation**: Update docstrings and README as you develop features

## Common Commands Summary

```bash
# Setup
pip install -e ".[dev]"

# Testing
pytest
pytest --cov=phylozoo --cov-report=html

# Examples
python examples/example_basic.py

# Development
mypy src/phylozoo
black src/phylozoo tests
ruff check src/phylozoo tests
```

