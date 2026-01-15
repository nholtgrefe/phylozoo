Advanced Topics
===============

This section covers advanced topics including custom exception handling, validation 
systems, and other advanced features of PhyloZoo.

Custom Exception Hierarchy
---------------------------

PhyloZoo uses a comprehensive custom exception hierarchy for clear, specific error 
messages. All exceptions inherit from ``PhyloZooError`` and are available at the 
package top-level:

.. code-block:: python

   from phylozoo import (
       PhyloZooError,           # Base exception
       PhyloZooValueError,      # Invalid values
       PhyloZooTypeError,       # Type errors
       PhyloZooNetworkError,    # Network-related errors
       PhyloZooIOError,         # I/O errors
       PhyloZooParseError,      # Parsing errors
       # ... and more
   )

Exception Categories
--------------------

**Value Errors**: ``PhyloZooValueError``
   Raised when a value is invalid (e.g., negative branch length, invalid gamma sum).

**Type Errors**: ``PhyloZooTypeError``
   Raised when a type is incorrect (e.g., passing a string where an integer is expected).

**Network Errors**: ``PhyloZooNetworkError`` and subclasses
   * ``PhyloZooNetworkStructureError``: Structural issues (cycles, connectivity)
   * ``PhyloZooNetworkDegreeError``: Degree constraint violations
   * ``PhyloZooNetworkAttributeError``: Attribute validation errors

**I/O Errors**: ``PhyloZooIOError`` and subclasses
   * ``PhyloZooParseError``: Parsing errors (malformed files)
   * ``PhyloZooFormatError``: Format-related errors

**Algorithm Errors**: ``PhyloZooAlgorithmError``
   Raised when algorithms encounter unsolvable problems.

**Visualization Errors**: ``PhyloZooVisualizationError`` and subclasses
   * ``PhyloZooLayoutError``: Layout computation errors
   * ``PhyloZooBackendError``: Backend-related errors
   * ``PhyloZooStateError``: Invalid visualization state

Example:

.. code-block:: python

   from phylozoo import PhyloZooNetworkStructureError
   
   try:
       network = DirectedPhyNetwork(edges=[(1, 2), (2, 1)])  # Cycle!
   except PhyloZooNetworkStructureError as e:
       print(f"Network structure error: {e}")

Validation System
-----------------

PhyloZoo includes a comprehensive validation system that ensures objects always 
represent valid phylogenetic structures.

Automatic Validation
--------------------

By default, network objects are validated upon creation:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   
   # This will raise PhyloZooNetworkStructureError due to cycle
   try:
       network = DirectedPhyNetwork(edges=[(1, 2), (2, 1)])
   except PhyloZooNetworkStructureError:
       print("Invalid network structure detected")

Validation Decorators
---------------------

The validation system uses decorators to control when validation occurs:

.. code-block:: python

   from phylozoo.utils.validation import validation_aware
   
   @validation_aware(allowed=["validate"], default=["validate"])
   class MyNetwork:
       def validate(self):
           # Validation logic
           pass

Disabling Validation
--------------------

For performance-critical code, validation can be temporarily disabled:

.. code-block:: python

   from phylozoo.utils.validation import no_validation
   
   with no_validation():
       # Operations that skip validation
       network = DirectedPhyNetwork(edges=[...])

**Warning**: Only disable validation when you are certain the operations produce 
valid structures. Invalid networks can cause errors in downstream operations.

Custom Validation
-----------------

You can add custom validation logic to network classes:

.. code-block:: python

   class MyNetwork(DirectedPhyNetwork):
       def _validate_custom_constraint(self):
           # Custom validation logic
           if some_condition:
               raise PhyloZooNetworkStructureError("Custom constraint violated")

Warnings
--------

PhyloZoo also provides custom warning classes:

.. code-block:: python

   from phylozoo import (
       PhyloZooWarning,
       PhyloZooEmptyNetworkWarning,
       PhyloZooSingleNodeNetworkWarning,
   )

Warnings are issued for non-critical issues (e.g., empty networks, single-node 
networks) that are technically valid but may not be useful for phylogenetic analysis.

Performance Considerations
---------------------------

**Caching**: PhyloZoo uses ``@cached_property`` for expensive computations. Properties 
are computed once and cached for subsequent access.

**NumPy**: Core data structures use NumPy arrays for efficient numerical operations.

**Numba**: Some computationally intensive algorithms use Numba JIT compilation for 
performance.

**Immutability**: Immutable objects enable safe caching and sharing without copying.

Best Practices
--------------

1. **Always handle exceptions**: Use try/except blocks for I/O and validation operations.

2. **Use type hints**: Type hints enable better IDE support and static type checking.

3. **Validate early**: Let PhyloZoo validate objects upon creation rather than 
   checking manually.

4. **Use appropriate exceptions**: Catch specific exception types rather than generic 
   ``Exception``.

5. **Disable validation carefully**: Only disable validation when absolutely necessary 
   for performance.
