I/O Operations
======================

This page describes the I/O protocol used in PhyloZoo and the basic operations
available on all I/O-capable classes.

IOMixin and I/O-capable classes
-------------------------------

All I/O in PhyloZoo goes through the :class:`~phylozoo.utils.io.mixin.IOMixin` protocol.
(If you are new to Python: a *mixin* is a class that provides methods to other
classes without being the main parent.) IOMixin adds ``save``, ``load``,
``to_string``, ``from_string``, ``convert``, and ``convert_string`` to any class that inherits from
it (see below). So instead of each class having its own I/O code, they all share the same
interface: you call ``network.save("file.enewick")`` or ``msa.load("file.fasta")``
in the same way regardless of the data type. PhyloZoo detects the format from the
file extension (e.g. ``.enewick``, ``.fasta``) and uses the registry to find the
right reader or writer for that class.

The following classes use IOMixin (directly or by inheritance):

- :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` — directed phylogenetic networks
- :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` — semi-directed phylogenetic networks
- :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` — directed multigraphs
- :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` — mixed directed/undirected multigraphs
- :class:`~phylozoo.core.split.splitsystem.SplitSystem` and :class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem` — split systems
- :class:`~phylozoo.core.distance.base.DistanceMatrix` — distance matrices
- :class:`~phylozoo.core.sequence.base.MSA` — multiple sequence alignments

Each class registers one or more formats with the central :class:`~phylozoo.utils.io.registry.FormatRegistry`; the
format(s) and extensions supported depend on the class. See :doc:`File formats
<formats/index>` for per-format details; see :doc:`Registering a new format <format_registry>` for how to register a new format for a class.

Basic operations
----------------

All I/O-capable classes support the following methods (defined on :class:`~phylozoo.utils.io.mixin.IOMixin`).
Methods may accept format-specific keyword arguments (e.g. ``triangle="UPPER"`` for NEXUS,
``root_location=...`` for semi-directed Newick). See the per-format pages for details.

Saving to a file
^^^^^^^^^^^^^^^^

The :meth:`~phylozoo.utils.io.mixin.IOMixin.save` method writes the object to a file. Format is detected from the file extension unless ``format`` is given. The ``overwrite`` parameter controls whether to overwrite the file if it already exists. 
Additional keyword arguments are passed to the format writer.

.. code-block:: python

   network.save("network.enewick")
   msa.save("alignment.fasta")
   dm.save("distances.nexus")

Loading from a file
^^^^^^^^^^^^^^^^^^^

The :meth:`~phylozoo.utils.io.mixin.IOMixin.load` method reads an object from a file. Format is detected from the extension unless ``format`` is given. If the extension is ambiguous or you use a non-standard extension, pass the format explicitly:

.. code-block:: python

   network = DirectedPhyNetwork.load("network.enewick")
   msa = MSA.load("alignment.fasta")
   dm = DistanceMatrix.load("distances.nexus")

   # Explicit format when extension is unclear
   network = DirectedPhyNetwork.load("network.txt", format="enewick")
   msa = MSA.load("data.txt", format="fasta")

Serializing to a string
^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~phylozoo.utils.io.mixin.IOMixin.to_string` method serializes the object to a string in the given format (default format for the class if not specified).

.. code-block:: python

   enewick_string = network.to_string(format="enewick")

Parsing from a string
^^^^^^^^^^^^^^^^^^^^^

The :meth:`~phylozoo.utils.io.mixin.IOMixin.from_string` method parses a string and returns an instance. Format must be specified or be the class default.

.. code-block:: python

   network = DirectedPhyNetwork.from_string(enewick_string, format="enewick")

Converting between file formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~phylozoo.utils.io.mixin.IOMixin.convert` method is a convenience wrapper that loads from one file and saves to another in one call. Formats are detected from extensions when not specified. Equivalent to ``cls.load(...).save(...)``; the object is created in memory during conversion.

.. code-block:: python

   MSA.convert("input.fasta", "output.nexus")

Converting between string formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :meth:`~phylozoo.utils.io.mixin.IOMixin.convert_string` method parses a string in one format and returns the serialized string in another format. Both formats must be specified. Equivalent to ``cls.from_string(content, format=input_format, **kwargs).to_string(format=output_format, **kwargs)``.

.. code-block:: python

   nexus_str = MSA.convert_string(fasta_content, input_format="fasta", output_format="nexus")

Error handling
--------------

I/O operations raise PhyloZoo exceptions when appropriate:

- :class:`~phylozoo.utils.exceptions.io.PhyloZooIOError` — General I/O errors (e.g. write failure).
- :class:`~phylozoo.utils.exceptions.io.PhyloZooParseError` — Parsing errors (malformed content).
- :class:`~phylozoo.utils.exceptions.io.PhyloZooFormatError` — Unsupported format, missing
  handler, or extension not recognized.
- :exc:`FileNotFoundError` — File not found (e.g. load from missing path).

.. code-block:: python

   from phylozoo.utils.exceptions import PhyloZooParseError

   try:
       network = DirectedPhyNetwork.load("malformed.enewick")
   except PhyloZooParseError as e:
       print(f"Parse error: {e}")

Best practices
--------------

- **Use standard file extensions** — PhyloZoo detects format from extensions. Use
  standard ones (e.g. ``.enewick``, ``.fasta``, ``.nexus``) when possible.

- **Specify format when needed** — If the extension is ambiguous or non-standard,
  pass the ``format`` parameter.

- **Handle exceptions** — Use try/except for I/O to handle parse errors and missing
  files gracefully.

- **Preserve metadata** — Some formats preserve more metadata (e.g. eNewick for
  branch lengths and gamma on networks). Choose a format that preserves what you need.

- **Use convert for convenience** — Use the ``convert`` class method to change
  format in one call instead of separate load and save.

See Also
--------

- :doc:`File formats <formats/index>` — Supported formats and format families
- :doc:`Registering a new format <format_registry>` — How to register a format
- :doc:`API Reference <../../../api/utils/io/index>` — The API reference for the I/O module