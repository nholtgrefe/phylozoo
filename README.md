# PhyloZoo

A phylogenetic analysis package for working with phylogenetic networks, trees, and related structures.

## Getting Started

### Installation

Install the package in development mode:

```bash
pip install -e .
```

This will install the package in editable mode, allowing you to make changes to the code without reinstalling.

### Development Installation

For development with all tools (testing, linting, type checking):

```bash
pip install -e ".[dev]"
```

Or just for testing:

```bash
pip install -e ".[test]"
```

## Usage

After installation, you can import the package:

```python
import phylozoo
```

## Testing

This project uses pytest for testing.

### Running Tests

To run all tests:

```bash
pytest
```

To run tests with coverage:

```bash
pytest --cov=phylozoo --cov-report=html
```

To run specific test files:

```bash
pytest tests/test_splits.py
```

To run tests matching a pattern:

```bash
pytest -k "test_split"
```

To run tests in verbose mode:

```bash
pytest -v
```

See `tests/README.md` for more testing information.


## Development

### Type Checking

This project uses type hints throughout. To check types:

```bash
mypy src/phylozoo
```

### Code Formatting

This project uses Black for code formatting:

```bash
black src/phylozoo tests
```

### Linting

This project uses Ruff for linting:

```bash
ruff check src/phylozoo tests
```

## License

MIT License

## Authors

- N. Holtgrefe (n.a.l.holtgrefe@tudelft.nl)

