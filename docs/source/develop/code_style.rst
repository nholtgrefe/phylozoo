Code Style Guide
================

This guide covers code style guidelines and conventions for PhyloZoo.

Type Hinting
------------

**All functions and classes must use type hints** for parameters, return types, and class 
attributes. Prefer modern built-in generics and PEP 604 unions (e.g., ``list[int]``, 
``dict[str, float]``, ``int | None``) now that we require Python 3.10+.

Example:

.. code-block:: python

   def process_network(
       network: DirectedPhyNetwork,
       threshold: float = 0.5,
   ) -> SemiDirectedPhyNetwork:
       """Process a network."""
       # Implementation
       pass

Use ``TypeVar`` for generic types when appropriate. Prefer explicit types over ``Any``—only 
use ``Any`` when absolutely necessary. Use ``-> None`` for functions that don't return values.

Docstrings
----------

**All public functions, classes, and methods must include NumPy-style docstrings** with 
detailed parameter descriptions, return values, exceptions, and examples.

Format:

.. code-block:: python

   def function_name(param1: Type, param2: Type) -> ReturnType:
       """
       Brief description of the function.
       
       Parameters
       ----------
       param1 : Type
           Description of param1.
       param2 : Type
           Description of param2.
       
       Returns
       -------
       ReturnType
           Description of what is returned.
       
       Raises
       ------
       ValueError
           When and why this exception is raised.
       
       Examples
       --------
       >>> example_usage()
       expected_output
       
       Notes
       -----
       Additional notes or implementation details.
       """

Include all sections that are relevant (Parameters, Returns, Raises, Examples, Notes, 
See Also, etc.). When including an example, also include it in the tests.

Naming Conventions
------------------

* **Classes**: Use ``PascalCase`` (e.g., ``DirectedPhyNetwork``, ``QuartetProfile``)
* **Functions and methods**: Use ``snake_case`` (e.g., ``compute_layout``, ``is_treechild``)
* **Constants**: Use ``UPPER_SNAKE_CASE`` (e.g., ``DEFAULT_RADIUS``)
* **Private/internal**: Use leading underscore (e.g., ``_internal_function``, ``_private_attr``)

Code Formatting
---------------

PhyloZoo uses **Black** for code formatting with a line length of 100 characters.

Format code before committing:

.. code-block:: bash

   black src/ tests/

Linting
-------

PhyloZoo uses **Ruff** for linting. Run linting checks:

.. code-block:: bash

   ruff check src/ tests/

Type Checking
-------------

PhyloZoo uses **mypy** for static type checking. Run type checks:

.. code-block:: bash

   mypy src/

Type checking configuration is in ``pyproject.toml``. The project uses Python 3.10+ with 
modern type hint features.

Code Quality Standards
----------------------

* **Immutability**: Core data structures should be immutable. Use frozen dataclasses or 
  similar patterns.
* **Validation**: Validate inputs early and provide clear error messages.
* **Error Handling**: Use the custom exception hierarchy from ``phylozoo.utils.exceptions``.
* **Performance**: Consider using NumPy for numerical operations and Numba for computationally 
  intensive algorithms when appropriate.
* **Documentation**: Keep docstrings up to date with code changes.

Example
-------

Here's an example of properly styled code:

.. code-block:: python

   from typing import TypeVar
   
   T = TypeVar('T')
   
   
   def compute_diversity(
       network: DirectedPhyNetwork,
       taxa: set[str],
       measure: DiversityMeasure,
   ) -> float:
       """
       Compute phylogenetic diversity for a set of taxa.
       
       Parameters
       ----------
       network : DirectedPhyNetwork
           The phylogenetic network.
       taxa : set[str]
           Set of taxon names to compute diversity for.
       measure : DiversityMeasure
           The diversity measure to use.
       
       Returns
       -------
       float
           The diversity value.
       
       Raises
       ------
       PhyloZooValueError
           If any taxon is not in the network.
       
       Examples
       --------
       >>> from phylozoo.panda import AllPathsDiversity
       >>> measure = AllPathsDiversity(network)
       >>> diversity = compute_diversity(network, {"A", "B", "C"}, measure)
       >>> diversity > 0
       True
       """
       if not taxa.issubset(network.taxa):
           raise PhyloZooValueError("All taxa must be in the network")
       
       return measure.compute_diversity(taxa)
