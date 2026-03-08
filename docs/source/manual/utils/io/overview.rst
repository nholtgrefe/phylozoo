Overview
========

PhyloZoo provides a unified I/O system for reading and writing phylogenetic data.
All I/O-capable classes use the same interface via :class:`~phylozoo.utils.io.IOMixin`
and a central :class:`~phylozoo.utils.io.FormatRegistry`.

This section is organized as follows:

- :doc:`I/O operations <operations>` — The IOMixin protocol, which classes
  support it, and the basic methods (save, load, to_string, from_string, convert).

- :doc:`File formats <formats/index>` — Supported formats and format families.
  Some formats form a *family*: the same file structure (e.g. NEXUS blocks) can
  represent different data types (distance matrices, alignments, split systems);
  the class you load or save with determines which variant is used.

- :doc:`Registering a new format <format_registry>` — How to register a new
  format for an existing or custom class using the FormatRegistry.
