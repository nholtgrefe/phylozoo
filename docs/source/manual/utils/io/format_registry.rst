Registering a New Format
========================

PhyloZoo's I/O system is driven by a central **FormatRegistry**. Each class that
supports I/O registers one or more formats with a reader and a writer. This page
explains how to register a new format for an existing class or for your own class.

Role of the Registry
--------------------

The :class:`phylozoo.utils.io.FormatRegistry` holds, for each (class, format name),
a **reader** (string → object) and a **writer** (object → string). When you call
``load``, ``save``, ``to_string``, or ``from_string`` on an I/O-capable class, the
registry is used to obtain the appropriate reader or writer. Format detection from
file extensions is also based on the registry (extensions are mapped to format names
per registration).

Registering a Format
--------------------

Use :meth:`FormatRegistry.register` to add a format for a class:

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

Parameters
^^^^^^^^^^

* **``obj_type``** — The class that supports this format (e.g. ``DirectedPhyNetwork``,
  ``MSA``).

* **``format``** — Format name string (e.g. ``'enewick'``, ``'nexus'``). This is the
  value users pass to ``format=`` in ``load``, ``save``, ``to_string``, ``from_string``.

* **``reader``** — Callable that parses a string and returns an instance of ``obj_type``.
  Signature: ``reader(string: str, **kwargs) -> obj_type``. Extra keyword arguments
  come from the caller (e.g. ``load(path, **kwargs)``).

* **``writer``** — Callable that serializes an instance to a string. Signature:
  ``writer(obj: obj_type, **kwargs) -> str``. Extra keyword arguments come from
  ``save`` or ``to_string``.

* **``extensions``** — Optional list of file extensions (e.g. ``['.nexus', '.nex']``).
  Each extension is associated with this format for format detection. Extensions are
  stored in lower case. If ``None``, the format cannot be detected from extension;
  users must pass ``format='myformat'`` explicitly.

* **``default``** — If ``True``, this format becomes the default for ``obj_type`` when
  no format is specified (e.g. ``to_string()`` with no ``format=`` argument). Only one
  default per class should be set.

Reader and Writer Contracts
---------------------------

* **Reader**: Receives the file content as a string (or the string passed to
  ``from_string``) and optional keyword arguments. Must return a single instance of
  the registered class. Should raise ``PhyloZooParseError`` (or similar) for
  malformed input.

* **Writer**: Receives an instance of the registered class and optional keyword
  arguments. Must return a string (which will be written to file or returned by
  ``to_string``). Should raise if the object cannot be serialized to this format.

Getting Handlers
----------------

You normally do not call these directly; ``IOMixin`` uses them internally. For
testing or custom code:

* **``FormatRegistry.get_reader(cls, format)``** — Returns the reader callable for
  that class and format. Raises ``PhyloZooFormatError`` if not registered.

* **``FormatRegistry.get_writer(cls, format)``** — Returns the writer callable.
  Raises ``PhyloZooFormatError`` if not registered.

* **``FormatRegistry.detect_format(filepath, obj_type)``** — Returns the format name
  detected from the file extension for the given class, or the class default if
  the extension is unknown. Raises ``PhyloZooFormatError`` if detection fails.

Example: Registering a Custom Format
-------------------------------------

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

.. seealso::
   * :doc:`overview` — Overview of the I/O section
   * :doc:`operations` — Using the I/O methods
   * :doc:`formats/index` — Built-in formats
   * :class:`phylozoo.utils.io.FormatRegistry` — API reference
