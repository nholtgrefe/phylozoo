Overview
========

The :mod:`phylozoo.utils` module provides supporting functionality for PhyloZoo:
exception handling, I/O operations, and validation. Import exceptions and
utilities from the submodules as needed.

.. code-block:: python

   from phylozoo.utils.exceptions import (
       PhyloZooError,
       PhyloZooValueError,
       PhyloZooNetworkError,
       PhyloZooParseError,
   )
   from phylozoo.utils.validation import no_validation, validation_aware

Submodules
----------

- :doc:`I/O <io/index>` — File I/O, format registry, and format-specific support
- :doc:`Exceptions <exceptions>` — Custom exception and warning hierarchy
- :doc:`Validation <validation>` — Network and generator validation

See Also
--------

- :doc:`API Reference <../../api/utils/index>` — Complete function and class reference
