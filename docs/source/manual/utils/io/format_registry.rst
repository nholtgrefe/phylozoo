Registering a New Format
========================

PhyloZoo's I/O system is driven by a central :class:`~phylozoo.utils.io.registry.FormatRegistry`.
Each class that supports I/O registers one or more formats with a reader and a writer.
This page explains how to register a new format for an existing class or for your own class.

Role of the Registry
--------------------

The :class:`~phylozoo.utils.io.registry.FormatRegistry` holds, for each (class, format name),
a **reader** (string → object) and a **writer** (object → string). When you call
:meth:`~phylozoo.utils.io.mixin.IOMixin.load`, :meth:`~phylozoo.utils.io.mixin.IOMixin.save`,
:meth:`~phylozoo.utils.io.mixin.IOMixin.to_string`, or :meth:`~phylozoo.utils.io.mixin.IOMixin.from_string`
on an I/O-capable class (i.e. a class that inherits from :class:`~phylozoo.utils.io.mixin.IOMixin`), the registry is used to obtain the appropriate reader or writer.
Format detection from file extensions is also based on the registry (extensions are mapped
to format names per registration).

Registering a Format
--------------------

Use :meth:`~phylozoo.utils.io.registry.FormatRegistry.register` to add a format for a class:

.. code-block:: python

   from phylozoo.utils.io import FormatRegistry

   FormatRegistry.register(
       MyClass,
       'myformat',
       reader=my_reader,
       writer=my_writer,
       extensions=['.myext', '.mef'],
       default=False,
   )

The :class:`~phylozoo.utils.io.registry.FormatRegistry` expects the following parameters:

**Reader**

The ``reader`` parameter is a callable that parses a string and returns an instance of the registered class.
Signature: ``reader(string: str, **kwargs) -> obj_type``.

**Writer**

The ``writer`` parameter is a callable that serializes an instance to a string. Signature:
``writer(obj: obj_type, **kwargs) -> str``.

**Extensions**

The ``extensions`` parameter is an optional list of file extensions (e.g. ``['.nexus', '.nex']``).
Each extension is associated with this format for format detection. Extensions are stored in lower
case. If ``None``, the format cannot be detected from extension; users must pass ``format='myformat'``
explicitly.

The ``default`` parameter, if ``True``, makes this format the default for the class when no format
is specified (e.g. :meth:`~phylozoo.utils.io.mixin.IOMixin.to_string` with no ``format=`` argument).
Only one default per class should be set.

Example
^^^^^^^

.. code-block:: python

   from phylozoo.utils.io import FormatRegistry
   from phylozoo.utils.exceptions import PhyloZooParseError

   def my_read(string: str, **kwargs):
       lines = string.strip().splitlines()
       if not lines:
           raise PhyloZooParseError("Empty input")
       # ... parse and return an instance of MyClass
       return MyClass(...)

   def my_write(obj: MyClass, **kwargs) -> str:
       # ... build string from obj
       return "..."

   FormatRegistry.register(
       MyClass, 'myformat',
       reader=my_read,
       writer=my_write,
       extensions=['.myext'],
       default=False,
   )

   # Then users can do:
   obj = MyClass.load("file.myext")  # format detected
   obj.save("out.myext")
   s = obj.to_string(format="myformat")

See Also
--------

- :doc:`General I/O operations <operations>` — Using the I/O methods
- :doc:`File formats <formats/index>` — Built-in formats
- :doc:`API Reference <../../../../api/utils/io/index>` — The API reference for the I/O module
