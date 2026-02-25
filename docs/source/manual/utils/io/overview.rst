I/O Documentation Overview
==========================

This section describes PhyloZoo's input/output system: how to use it, how to extend it,
and the file formats that are supported.

Pages in this section
---------------------

* **:doc:`operations`** — How to use the I/O methods: ``save``, ``load``, ``to_string``,
  ``from_string``, and ``convert``. Format detection, error handling, and best practices.

* **:doc:`format_registry`** — How to register a new file format for a class using
  ``FormatRegistry`` (reader, writer, extensions, default format).

* **:doc:`formats/index`** — One page per format (or format family), with structure,
  which classes use it, file extensions, and examples:

  * :doc:`NEXUS <formats/nexus>` — Block-structured format (distance matrices, MSAs, split systems)
  * :doc:`PHYLIP <formats/phylip>` — Distance matrices
  * :doc:`Newick / eNewick <formats/newick>` — Trees and directed networks
  * :doc:`DOT <formats/dot>` — Graphviz and PhyloZoo-DOT for graphs and networks
  * :doc:`FASTA <formats/fasta>` — Sequence alignments
  * :doc:`CSV <formats/csv>` — Distance matrices (tabular)
  * :doc:`Edge list <formats/edgelist>` — Simple graph edge lists
