Exceptions
==========

The :mod:`phylozoo.utils.exceptions` module provides a custom exception hierarchy for clear, specific error and warningmessages.

All exceptions and functions on this page can be imported from the exceptions module:

.. code-block:: python

   from phylozoo.utils.exceptions import (
       PhyloZooError,
       PhyloZooValueError,
       PhyloZooNetworkError,
       PhyloZooParseError,
       PhyloZooWarning,
       warn_on_keyword,
       warn_on_none_value,
   )

Errors
------

Errors are raised when something goes wrong in PhyloZoo.


Base Error Class
^^^^^^^^^^^^^^^^^

:class:`~phylozoo.utils.exceptions.base.PhyloZooError` is the base exception for all PhyloZoo errors.
All custom exceptions inherit from this class. Catch it to handle any PhyloZoo-specific error.

.. code-block:: python

   try:
       network = DirectedPhyNetwork.load("network.enewick")
   except PhyloZooError as e:
       print(f"PhyloZoo error: {e}")

General-Purpose Errors
^^^^^^^^^^^^^^^^^^^^^^^

General-purpose errors inherit from both
:class:`~phylozoo.utils.exceptions.base.PhyloZooError` and the corresponding Python built-in,
so code that catches, e.g., :exc:`ValueError` or :exc:`TypeError` still works:

- :class:`~phylozoo.utils.exceptions.general.PhyloZooValueError` — also inherits from :exc:`ValueError`
- :class:`~phylozoo.utils.exceptions.general.PhyloZooTypeError` — also inherits from :exc:`TypeError`
- :class:`~phylozoo.utils.exceptions.general.PhyloZooRuntimeError` — also inherits from :exc:`RuntimeError`
- :class:`~phylozoo.utils.exceptions.general.PhyloZooNotImplementedError` — also inherits from :exc:`NotImplementedError`
- :class:`~phylozoo.utils.exceptions.general.PhyloZooImportError` — also inherits from :exc:`ImportError`
- :class:`~phylozoo.utils.exceptions.general.PhyloZooAttributeError` — also inherits from :exc:`AttributeError`

Domain-specific Errors
^^^^^^^^^^^^^^^^^^^^^^^

Each domain has a base exception; subclasses provide finer-grained handling. See the
:doc:`Exceptions API <../../api/utils/exceptions>` for all subclasses.

- :class:`~phylozoo.utils.exceptions.network.PhyloZooNetworkError` — network structure, degree, and attribute errors
- :class:`~phylozoo.utils.exceptions.io.PhyloZooIOError` — I/O, parse, and format errors
- :class:`~phylozoo.utils.exceptions.algorithm.PhyloZooAlgorithmError` — algorithm-related errors
- :class:`~phylozoo.utils.exceptions.visualization.PhyloZooVisualizationError` — layout, backend, and state errors
- :class:`~phylozoo.utils.exceptions.generator.PhyloZooGeneratorError` — network generator structure and degree errors

Warnings
--------

Warnings are raised when the user does something that is not strictly wrong, but might cause unexpected behavior or issues.

Base Warning Class
^^^^^^^^^^^^^^^^^^

:class:`~phylozoo.utils.exceptions.base.PhyloZooWarning` is the base exception for all PhyloZoo warnings.
All custom warnings inherit from this class. Catch it to handle any PhyloZoo-specific warning.

.. code-block:: python

   try:
       network = DirectedPhyNetwork.load("network.enewick")
   except PhyloZooWarning as e:
       print(f"PhyloZoo warning: {e}")

Specific Warnings
^^^^^^^^^^^^^^^^^

- :class:`~phylozoo.utils.exceptions.warning.PhyloZooIdentifierWarning` — warn if a value is a Python keyword (e.g. used as an identifier)
- :class:`~phylozoo.utils.exceptions.warning.PhyloZooEmptyNetworkWarning` — warn if a network is empty
- :class:`~phylozoo.utils.exceptions.warning.PhyloZooSingleNodeNetworkWarning` — warn if a network has only one node

Warning Utilities
^^^^^^^^^^^^^^^^^

The exceptions module provides helpers that emit PhyloZoo warnings:

- :func:`~phylozoo.utils.exceptions.utils.warn_on_keyword` — warn if a value is a Python keyword (e.g. used as an identifier)
- :func:`~phylozoo.utils.exceptions.utils.warn_on_none_value` — warn if a value is ``None`` (Python keyword)

.. code-block:: python

   from phylozoo.utils.exceptions import warn_on_keyword, warn_on_none_value

   warn_on_keyword("for", "Identifier")
   warn_on_none_value(None, "Attribute 'weight'")


See Also
--------

- :doc:`Exceptions API <../../api/utils/exceptions>` — full list of exception and warning classes.
- :doc:`I/O <io/index>` — I/O operations and parse/format errors.
- :doc:`Validation <validation>` — network validation.
