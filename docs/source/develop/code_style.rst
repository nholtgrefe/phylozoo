Code Style Guide
================

This guide covers code style guidelines and conventions for PhyloZoo.

Code Conventions
----------------

Type Hinting
~~~~~~~~~~~~

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
~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~

* **Classes**: Use ``PascalCase`` (e.g., ``DirectedPhyNetwork``, ``QuartetProfile``)
* **Functions and methods**: Use ``snake_case`` (e.g., ``compute_layout``, ``is_treechild``)
* **Constants**: Use ``UPPER_SNAKE_CASE`` (e.g., ``DEFAULT_RADIUS``)
* **Private/internal**: Use leading underscore (e.g., ``_internal_function``, ``_private_attr``)

Development Tools
-----------------

PhyloZoo relies on standard tools for code quality checks and automation: Black (formatting),
Ruff (linting), and mypy (type checking).

Installation
~~~~~~~~~~~~

Install the tools directly:

.. code-block:: bash

   pip install black ruff mypy

or via the development extra:

.. code-block:: bash

   pip install -e ".[dev]"

Code Formatting (Black)
~~~~~~~~~~~~~~~~~~~~~~~

PhyloZoo uses **Black** for code formatting with a line length of 100 characters.

Format code before committing:

.. code-block:: bash

   black src/ tests/

Linting (Ruff)
~~~~~~~~~~~~~~

PhyloZoo uses **Ruff** for linting. Run linting checks:

.. code-block:: bash

   ruff check src/ tests/

Type Checking (mypy)
~~~~~~~~~~~~~~~~~~~~

PhyloZoo uses **mypy** for static type checking. Run type checks:

.. code-block:: bash

   mypy src/

Type checking configuration is in ``pyproject.toml``. The project uses Python 3.10+ with 
modern type hint features.

Example
-------

Here's an example of properly styled code:

.. code-block:: python

   from typing import TypeVar
   from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
   
   T = TypeVar('T')
   
   
   def induced_subnetwork(
       network: SemiDirectedPhyNetwork,
       taxa: list[str],
   ) -> SemiDirectedPhyNetwork:
       """
       Extract the subnetwork induced by a subset of taxa.
       
       Parameters
       ----------
       network : SemiDirectedPhyNetwork
           The phylogenetic network.
       taxa : list[str]
           Subset of taxon labels to induce the subnetwork on.
       
       Returns
       -------
       SemiDirectedPhyNetwork
           The induced subnetwork.
       
       Raises
       ------
       PhyloZooValueError
           If any taxon is not in the network.
       
       Examples
       --------
       >>> from phylozoo.core.network.sdnetwork import derivations
       >>> network = ...  # SemiDirectedPhyNetwork
       >>> sub = derivations.subnetwork(network, taxa=["A", "B", "C"])
       >>> len(sub.taxa) <= 3
       True
       """
       from phylozoo.core.network.sdnetwork import derivations
       return derivations.subnetwork(network, taxa)
