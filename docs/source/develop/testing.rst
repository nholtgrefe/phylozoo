Testing Guide
==============

This guide covers testing practices and guidelines for PhyloZoo.

Setup
-----

Testing Dependencies
~~~~~~~~~~~~~~~~~~~~

PhyloZoo uses a small set of testing tools:

* `pytest <https://docs.pytest.org/>`_ >= 7.0.0
* `pytest-cov <https://pytest-cov.readthedocs.io/>`_ >= 4.0.0

You can install them directly:

.. code-block:: bash

   pip install pytest pytest-cov

or as part of the development extra when working on PhyloZoo itself:

.. code-block:: bash

   pip install -e ".[dev]"

Test Organization
-----------------

Test Framework
~~~~~~~~~~~~~~

PhyloZoo uses `pytest` as the testing framework. Tests are located in the ``tests/`` directory
and mirror the source code structure.

Test Structure
~~~~~~~~~~~~~~

Tests are organized to mirror the source code structure:

* `tests/core/` - Tests for core modules (networks, quartets, splits, sequences, distance)
* `tests/viz/` - Tests for visualization
* `tests/utils/` - Tests for utility modules
* `tests/conftest.py` - Shared fixtures and pytest configuration

Running Tests
-------------

**Run all tests:**

.. code-block:: bash

   pytest

**Run tests with coverage:**

.. code-block:: bash

   pytest --cov=phylozoo --cov-report=html

**Run specific test files:**

.. code-block:: bash

   pytest tests/core/network/test_dnetwork.py

**Run tests matching a pattern:**

.. code-block:: bash

   pytest -k "test_split"

**Run tests and show print statements:**

.. code-block:: bash

   pytest -s

**Run tests in verbose mode:**

.. code-block:: bash

   pytest -v

Writing Tests
-------------

Guidelines
~~~~~~~~~~

When adding new tests:

1. **Create a test file** named `test_<module_name>.py` in the appropriate subdirectory
2. **Import the module** you're testing
3. **Use descriptive test function names** starting with `test_`
4. **Use pytest fixtures** from `conftest.py` when appropriate
5. **Add docstrings** to test classes and functions explaining what they test
6. **Use type hints** for test functions

Example:

.. code-block:: python

   def test_my_function() -> None:
       """Test that my_function works correctly."""
       result = my_function(input)
       assert result == expected_output

Test Classes
~~~~~~~~~~~~

Group related tests in test classes:

.. code-block:: python

   class TestMyClass:
       """Test cases for MyClass."""
       
       def test_basic_functionality(self) -> None:
           """Test basic functionality."""
           # Test code here
           pass
       
       def test_edge_case(self) -> None:
           """Test edge case handling."""
           # Test code here
           pass

Fixtures
~~~~~~~~

Shared fixtures are defined in `tests/conftest.py`. Use fixtures for common test data:

.. code-block:: python

   def test_with_fixture(sample_network) -> None:
       """Test using a fixture."""
       assert sample_network.num_nodes > 0

Best Practices
~~~~~~~~~~~~~~

* **Test both success and failure cases**: Test that functions work correctly and handle 
  errors appropriately
* **Use descriptive assertions**: Include clear error messages in assertions when possible
* **Test edge cases**: Test with empty inputs, single elements, boundary conditions
* **Keep tests independent**: Each test should be able to run independently
* **Use fixtures for common setup**: Avoid duplicating test setup code
