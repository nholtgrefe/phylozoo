Exceptions
==========

PhyloZoo uses a comprehensive custom exception hierarchy for clear, specific error messages. 
All exceptions inherit from ``PhyloZooError`` and are available at the package top-level, 
allowing you to catch PhyloZoo-specific errors easily.

Exception Hierarchy
-------------------

All PhyloZoo exceptions inherit from ``PhyloZooError``:

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

You can catch all PhyloZoo errors with a single except clause:

.. code-block:: python

   try:
       # Some PhyloZoo operation
       network = DirectedPhyNetwork.load("network.enewick")
   except PhyloZooError as e:
       print(f"PhyloZoo error: {e}")

General Purpose Exceptions
---------------------------

**PhyloZooError**
   Base exception for all PhyloZoo errors. All custom exceptions inherit from this class.

**PhyloZooValueError**
   Raised when a value is invalid (e.g., negative branch length, invalid gamma sum, 
   out-of-range index). Inherits from both ``PhyloZooError`` and built-in ``ValueError``.

**PhyloZooTypeError**
   Raised when a type is incorrect (e.g., passing a string where an integer is expected). 
   Inherits from both ``PhyloZooError`` and built-in ``TypeError``.

**PhyloZooRuntimeError**
   Raised when a runtime error occurs. Inherits from both ``PhyloZooError`` and built-in 
   ``RuntimeError``.

**PhyloZooNotImplementedError**
   Raised when a feature is not yet implemented or not available for the given input. 
   Inherits from both ``PhyloZooError`` and built-in ``NotImplementedError``.

**PhyloZooAttributeError**
   Raised when an attribute access fails. Inherits from both ``PhyloZooError`` and built-in 
   ``AttributeError``.

Network Domain Exceptions
-------------------------

**PhyloZooNetworkError**
   Base exception for network-related errors. Inherits from ``PhyloZooValueError``. All 
   network-specific errors inherit from this class.

**PhyloZooNetworkStructureError**
   Raised when network structure is invalid:
   * Directed cycles in directed networks
   * Undirected cycles in semi-directed networks
   * Disconnected networks
   * Multiple root nodes
   * Invalid node degrees

**PhyloZooNetworkDegreeError**
   Raised when node degree constraints are violated:
   * Internal nodes with insufficient degree
   * Hybrid nodes with invalid in-degree/out-degree
   * Tree nodes with invalid in-degree/out-degree

**PhyloZooNetworkAttributeError**
   Raised when network attribute validation fails:
   * Invalid branch lengths
   * Invalid gamma values (not summing to 1.0 for hybrid nodes)
   * Invalid bootstrap values (not in [0.0, 1.0])
   * Parallel edges with different branch lengths

I/O Domain Exceptions
---------------------

**PhyloZooIOError**
   Base exception for I/O-related errors. Inherits from both ``PhyloZooError`` and built-in 
   ``IOError``. Raised for general I/O problems.

**PhyloZooParseError**
   Raised when parsing fails (malformed files, invalid format syntax). Inherits from 
   ``PhyloZooIOError``. Used for:
   * eNewick parsing errors
   * Newick parsing errors
   * NEXUS parsing errors
   * Other format parsing errors

**PhyloZooFormatError**
   Raised when format-related errors occur (unsupported format, missing format handler, 
   format mismatch). Inherits from ``PhyloZooIOError``.

Algorithm Domain Exceptions
---------------------------

**PhyloZooAlgorithmError**
   Base exception for algorithm-related errors. Raised when algorithms encounter unsolvable 
   problems or invalid input configurations.

Visualization Domain Exceptions
--------------------------------

**PhyloZooVisualizationError**
   Base exception for visualization-related errors.

**PhyloZooLayoutError**
   Raised when layout computation fails (e.g., non-planar graph for planar layout, 
   invalid layout parameters).

**PhyloZooBackendError**
   Raised when visualization backend errors occur.

**PhyloZooStateError**
   Raised when visualization state is invalid.

Common Exception Patterns
-------------------------

**Catching Specific Exceptions**

.. code-block:: python

   from phylozoo import (
       PhyloZooNetworkStructureError,
       PhyloZooParseError,
       PhyloZooValueError
   )
   
   try:
       network = DirectedPhyNetwork(edges=[(1, 2), (2, 1)])  # Cycle!
   except PhyloZooNetworkStructureError as e:
       print(f"Network structure error: {e}")
   except PhyloZooValueError as e:
       print(f"Value error: {e}")

**Catching All PhyloZoo Errors**

.. code-block:: python

   from phylozoo import PhyloZooError
   
   try:
       network = DirectedPhyNetwork.load("network.enewick")
       result = some_operation(network)
   except PhyloZooError as e:
       print(f"PhyloZoo error occurred: {e}")
       # Handle error appropriately

**Catching I/O Errors**

.. code-block:: python

   from phylozoo import PhyloZooIOError, PhyloZooParseError
   
   try:
       network = DirectedPhyNetwork.load("malformed.enewick")
   except PhyloZooParseError as e:
       print(f"Parse error: {e}")
   except PhyloZooIOError as e:
       print(f"I/O error: {e}")

**Catching Network Errors**

.. code-block:: python

   from phylozoo import PhyloZooNetworkError
   
   try:
       network = DirectedPhyNetwork(edges=[...])
   except PhyloZooNetworkError as e:
       print(f"Network error: {e}")
       # All network-related errors (structure, degree, attribute) are caught

Best Practices
--------------

1. **Catch specific exceptions**: Catch the most specific exception type possible. This 
   allows you to handle different error types appropriately.

2. **Use PhyloZooError for general handling**: When you want to catch all PhyloZoo errors 
   but don't need specific handling, catch ``PhyloZooError``.

3. **Preserve exception information**: Always include the exception message when logging 
   or displaying errors. PhyloZoo exceptions include detailed error messages.

4. **Handle expected errors**: Some operations may legitimately raise exceptions (e.g., 
   ``PhyloZooNotImplementedError`` when exact solution is not available). Handle these 
   gracefully.

5. **Don't catch generic Exception**: Avoid catching generic ``Exception`` unless necessary. 
   Use ``PhyloZooError`` or more specific exceptions instead.

6. **Check exception types**: Use ``isinstance()`` to check exception types if you need 
   to handle multiple exception types in a single except clause.

Example: Comprehensive Error Handling
--------------------------------------

.. code-block:: python

   from phylozoo import (
       PhyloZooError,
       PhyloZooNetworkStructureError,
       PhyloZooParseError,
       PhyloZooNotImplementedError
   )
   from phylozoo import DirectedPhyNetwork
   
   try:
       # Load network
       network = DirectedPhyNetwork.load("network.enewick")
       # ... use network ...
   except PhyloZooParseError as e:
       print(f"Failed to parse network file: {e}")
   except PhyloZooNetworkStructureError as e:
       print(f"Invalid network structure: {e}")
   except PhyloZooError as e:
       print(f"PhyloZoo error: {e}")
   except Exception as e:
       print(f"Unexpected error: {e}")

.. seealso::
   For more information on:
   * I/O operations: :doc:`I/O <io>`
   * Network validation: :doc:`Validation <advanced/validation>`
   * Error handling in specific modules: See module-specific documentation
