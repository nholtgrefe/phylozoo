I/O Operations
=============

This page describes the I/O methods provided by PhyloZoo: saving and loading files,
converting to and from strings, and converting between formats. All I/O-capable classes
use the same interface via ``IOMixin``.

I/O System
---------

All I/O-capable classes in PhyloZoo use the ``IOMixin`` protocol, which provides a
consistent interface for loading, saving, and converting between formats. The system
automatically detects file formats based on file extensions and supports explicit
format specification.

Basic Operations
----------------

All I/O-capable classes support these methods:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork, MSA, DistanceMatrix

   # Save to file (auto-detects format from extension)
   network.save("network.enewick")
   msa.save("alignment.fasta")
   dm.save("distances.nexus")

   # Load from file (auto-detects format)
   network = DirectedPhyNetwork.load("network.enewick")
   msa = MSA.load("alignment.fasta")
   dm = DistanceMatrix.load("distances.nexus")

   # Convert to string
   enewick_string = network.to_string(format="enewick")

   # Parse from string
   network = DirectedPhyNetwork.from_string(enewick_string, format="enewick")

   # Convert between formats (file to file, without loading into memory)
   MSA.convert("input.fasta", "output.nexus")

Method summary
^^^^^^^^^^^^^^

* **``save(filepath, format=None, overwrite=True, **kwargs)``** — Write the object to a
  file. Format is detected from the file extension unless ``format`` is given.
  Additional keyword arguments are passed to the format writer.

* **``load(filepath, format=None, **kwargs)``** — Class method. Read an object from a
  file. Format is detected from the extension unless ``format`` is given.

* **``to_string(format=None, **kwargs)``** — Serialize the object to a string in the
  given format (default format for the class if not specified).

* **``from_string(string, format=None, **kwargs)``** — Class method. Parse a string
  and return an instance. Format must be specified or be the class default.

* **``convert(input_path, output_path, input_format=None, output_format=None, **kwargs)``**
  — Class method. Read from one file and write to another. Formats are detected from
  extensions when not specified. Useful for large files when you do not need the
  object in memory.

Format Detection
----------------

PhyloZoo automatically detects file formats based on file extensions. If the extension
is ambiguous or you use a non-standard extension, specify the format explicitly:

.. code-block:: python

   # Explicit format specification
   network = DirectedPhyNetwork.load("network.txt", format="enewick")
   msa = MSA.load("data.txt", format="fasta")

Error Handling
--------------

I/O operations raise appropriate exceptions when files cannot be read or written:

* **``PhyloZooIOError``** — General I/O errors (e.g. write failure).

* **``PhyloZooParseError``** — Parsing errors (malformed file content).

* **``PhyloZooFormatError``** — Format-related errors (unsupported format, missing
  format handler, or extension not recognized).

* **``FileNotFoundError``** — File not found (e.g. load from missing path).

Example:

.. code-block:: python

   from phylozoo.utils.exceptions import PhyloZooParseError

   try:
       network = DirectedPhyNetwork.load("malformed.enewick")
   except PhyloZooParseError as e:
       print(f"Parse error: {e}")

Best Practices
--------------

1. **Use appropriate file extensions** — PhyloZoo auto-detects formats from extensions.
   Use standard extensions (e.g. ``.enewick``, ``.fasta``, ``.nexus``) for best results.

2. **Specify format explicitly when needed** — If the extension is ambiguous or
   non-standard, pass the ``format`` parameter.

3. **Handle exceptions** — Use try/except for I/O operations to handle parsing errors
   and missing files gracefully.

4. **Preserve metadata** — Some formats preserve more metadata than others. Choose a
   format that preserves the attributes you need (e.g. eNewick for branch lengths and
   gamma values on networks).

5. **Convert between formats** — Use the ``convert`` class method to convert between
   formats without loading into memory when working with large files.

.. seealso::
   * :doc:`overview` — Overview of this I/O section
   * :doc:`format_registry` — Registering a new format
   * :doc:`formats/index` — Supported file formats
   * :doc:`../exceptions` — Exception types
