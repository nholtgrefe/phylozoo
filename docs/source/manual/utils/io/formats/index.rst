File Formats
============

This page lists all currently supported file formats in PhyloZoo. Each format is
documented on its own page with structure, which classes use it, file extensions,
and examples.

Some formats act as *families*: the same file structure or naming can represent
different data types; the class you use for load/save determines which variant
is read or written. Others are single-purpose.

- :doc:`NEXUS <nexus>` — Block-structured (``BEGIN ... END;``). Different block
  types hold distance matrices, character matrices (MSA), or split systems.

- :doc:`PHYLIP <phylip>` — Matrix-style layout (first line = size, then labeled
  rows). Used for distance matrices.

- :doc:`eNewick <enewick>` — Standard format for trees and networks.
  Used by both directed and semi-directed networks; eNewick extends Newick with
  hybrid nodes and edge attributes.

- :doc:`DOT <dot>` — Graphviz format for directed graphs. PhyloZoo-DOT (``pzdot``)
  extends it for semi-directed networks and mixed graphs.

- :doc:`FASTA <fasta>` — Simple text format for sequence alignments (header
  ``>``, then sequence lines).

- :doc:`CSV <csv>` — Comma-separated values for distance matrices, with optional
  header row and custom delimiters.

- :doc:`Edge list <edgelist>` — One edge per line (two node IDs). For directed
  multigraphs.

.. note::
   I/O for each class lives in the specific module (e.g. ``core/distance/io`` for
   :class:`~phylozoo.core.distance.base.DistanceMatrix`, ``core/sequence/io`` for
   :class:`~phylozoo.core.sequence.base.MSA`). See the :doc:`API Reference <../../../../api/utils/io/index>` for the full I/O
   module layout.

   Format families may share helper functions for parsing or writing (e.g. NEXUS blocks, PHYLIP matrix layout). These shared
   utilities live in :mod:`phylozoo.utils.io.format_utils` and are used by the
   class-specific I/O modules. See :doc:`Format utilities <../../../../api/utils/io/format_utils/index>`
   in the API reference.

.. toctree::
   :hidden:
   :maxdepth: 1

   nexus
   phylip
   enewick
   dot
   fasta
   csv
   edgelist
