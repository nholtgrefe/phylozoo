# Tests

This directory contains tests for the phylozoo package.

## Running Tests

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

To run tests and show print statements:

```bash
pytest -s
```

To run tests in verbose mode:

```bash
pytest -v
```

To run only fast tests (excluding slow markers):

```bash
pytest -m "not slow"
```

## Test Structure

Tests are organized to mirror the source code structure:
- `test_splits.py` - Tests for the splits module
- `test_utils.py` - Tests for utility modules
- `test_networks.py` - Tests for network modules
- `conftest.py` - Shared fixtures and pytest configuration

## Writing New Tests

When adding new tests:

1. Create a test file named `test_<module_name>.py`
2. Import the module you're testing
3. Use descriptive test function names starting with `test_`
4. Use pytest fixtures from `conftest.py` when appropriate
5. Add docstrings to test classes and functions
6. Use type hints for test functions

Example:

```python
def test_my_function() -> None:
    """Test that my_function works correctly."""
    result = my_function(input)
    assert result == expected_output
```

## Test Markers

The following markers are available:
- `@pytest.mark.slow` - Marks tests as slow (deselect with '-m "not slow"')
- `@pytest.mark.integration` - Marks tests as integration tests

