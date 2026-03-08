Validation
==========

The :mod:`phylozoo.utils.validation` module contains the utilities used to validate PhyloZoo objects and to disable validation for performance-critical code.

All functions and classes on this page can be imported from the validation module:

.. code-block:: python

   from phylozoo.utils.validation import no_validation, validation_aware

What is validation?
-------------------

PhyloZoo performs validation on objects to ensure they satisfy required structural and attribute rules. 
For example, these checks cover connectivity, the absence of forbidden cycles, correct node degrees (such as having exactly one root and leaves with in-degree 1), and valid edge attributes (for example, ensuring hybrid-edge γ values lie in [0, 1] and sum to 1.0 at hybrid nodes).

By default, validation runs automatically at the end of object construction so that invalid objects are rejected as soon as they are created. 
For some classes, validation is expensive, and it can therefore be disabled; for performance, or when building objects in stages.

If a check fails, a domain‑specific exception is raised (such as :class:`~phylozoo.utils.exceptions.network.PhyloZooNetworkStructureError` or :class:`~phylozoo.utils.exceptions.network.PhyloZooNetworkDegreeError`). 
See :doc:`Exceptions <exceptions>` for the full exception hierarchy.
This page explains which classes use the validation system, how to enable or disable it, and how to make your own classes validation‑aware.

Which classes are validation-aware?
--------------------------------

The following classes are validation-aware and run a validation method that is recognized by the validation system:

- :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` — directed phylogenetic networks
- :class:`~phylozoo.core.network.sdnetwork.base.MixedPhyNetwork` — mixed directed/undirected networks (base for semi-directed)
- :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` — semi-directed phylogenetic networks
- :class:`~phylozoo.core.network.dnetwork.generator.base.DirectedGenerator` — directed level-*k* generators
- :class:`~phylozoo.core.network.sdnetwork.generator.base.SemiDirectedGenerator` — semi-directed level-*k* generators

For the exact checks each class performs, see their documentation in the :doc:`API reference <../../api/index>`.

Note that some other classes also validate input but they are are not validation-aware.
This design choice was made for classes where validation is inexpensive.

Disabling validation
--------------------

Use the :func:`~phylozoo.utils.validation.no_validation` context manager to
temporarily disable validation. Inside the block, calls to :meth:`validate` (and
any other methods you configured) are skipped for the matching classes and
methods, so construction no longer runs those checks.

**Disable all validation (default methods for all validation-aware classes)**

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.utils.validation import no_validation

   with no_validation():
       net = DirectedPhyNetwork(edges=[(1, 2), (2, 3)])  # validate() not run

**Disable only for specific classes**

Use the ``classes`` argument with class names or ``fnmatch`` patterns. Only
those classes skip validation; others still run it.

.. code-block:: python

   from phylozoo.core.network.dnetwork.generator.base import DirectedGenerator
   from phylozoo.utils.validation import no_validation

   with no_validation(classes=["DirectedGenerator"]):
       gen = DirectedGenerator(...)        # validate() suppressed
       net = DirectedPhyNetwork(...)       # validate() still runs

   with no_validation(classes=["*Generator"]):  # any class whose name ends with "Generator"
       gen = DirectedGenerator(...)         # validate() suppressed

**Disable only specific methods**

Use the ``methods`` argument with method names or `fnmatch` patterns. The
class’s ``default`` list (see below) is ignored when you pass ``methods``.

.. code-block:: python

   from phylozoo.utils.validation import no_validation

   with no_validation(methods=["validate"]):
       net = DirectedPhyNetwork(...)  # only validate() is suppressed

   with no_validation(methods=["validate", "_validate_*"]):
       net = DirectedPhyNetwork(...)  # validate() and all _validate_* helpers suppressed

**Combine class and method filters**

.. code-block:: python

   with no_validation(classes=["*Generator"], methods=["validate"]):
       gen = DirectedGenerator(...)  # validate() suppressed only for generator classes

.. note::
    Nested blocks stack: inner ``no_validation`` settings apply in addition to
    outer ones. After exiting the context, validation runs normally again.

.. warning::
   Only turn off validation when you are sure the objects you build are valid.
   Invalid objects can cause confusing errors in later code.


Example: Batch creation of networks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When building many networks with a particular structures, which are known or tested to be valid, it can be beneficial to disable validation for performance.
After constructing the networks, you can still validate them individually if needed.

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.utils.validation import no_validation
   from phylozoo.utils.exceptions import PhyloZooNetworkError

   with no_validation():
       networks = [
           DirectedPhyNetwork(edges=[(i, i + 1)])
           for i in range(100)
       ]

   for net in networks:
       try:
           net.validate()
       except PhyloZooNetworkError as e:
           print(f"Invalid network: {net}, error: {e}")


Using the decorator for your own classes
----------------------------------------

To make a class participate in this system, decorate it with
:func:`~phylozoo.utils.validation.validation_aware`. You choose which methods
can be suppressed and which are suppressed by default when
:func:`~phylozoo.utils.validation.no_validation` is used with no
``methods`` argument.

**Parameters**

- **allowed** — Method names or ``fnmatch`` patterns that may be suppressed
  (e.g. ``["validate", "_validate_*"]``). Only these methods are wrapped.
- **default** — Subset of ``allowed`` that is suppressed when
  :func:`~phylozoo.utils.validation.no_validation` is used without
  ``methods``. Must match at least one ``allowed`` pattern or
  :class:`~phylozoo.utils.exceptions.general.PhyloZooValueError` is raised.

**Example**

.. code-block:: python

   from phylozoo.utils.validation import validation_aware, no_validation

   @validation_aware(allowed=["validate", "_validate_*"], default=["validate"])
   class MyNetwork:
       def __init__(self):
           self.validate()  # runs unless inside no_validation()

       def validate(self):
           self._validate_structure()

       def _validate_structure(self):
           # your checks; raise PhyloZoo*Error on failure
           pass

   # Normal use: validation runs
   obj = MyNetwork()

   # Suppress only validate (default for this class)
   with no_validation():
       obj2 = MyNetwork()  # validate() and _validate_* not run

   # Suppress only a specific method for all classes
   with no_validation(methods=["_validate_structure"]):
       obj.validate()  # validate() runs, _validate_structure() is skipped


See Also
--------

- :doc:`Exceptions <exceptions>` — Exception and warning classes raised during validation.
- :doc:`Directed Network Class <../core/networks/directed/directed_network_class>` — DirectedPhyNetwork construction and validation.
- :doc:`Semi-Directed Network Class <../core/networks/semi_directed/semi_directed_network_class>` — SemiDirectedPhyNetwork and MixedPhyNetwork construction and validation.
- :doc:`Directed Generator <../core/networks/directed/directed_generator>` — DirectedGenerator construction and validation.
- :doc:`Semi-Directed Generator <../core/networks/semi_directed/generators>` — SemiDirectedGenerator construction and validation.
- :doc:`Validation API <../../api/utils/index>`
