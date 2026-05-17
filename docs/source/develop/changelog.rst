Changelog
=========

Version History
---------------

Version 0.2.0
~~~~~~~~~~~~~

Added
^^^^^

* New ``phylozoo.core.triplet`` module providing :class:`~phylozoo.core.triplet.base.Triplet`,
  :class:`~phylozoo.core.triplet.tprofile.TripletProfile`, and
  :class:`~phylozoo.core.triplet.tprofileset.TripletProfileSet` for working with rooted
  three-taxon trees (characterised by a trivial 1|2 split ``a|bc`` or by an unresolved
  3-taxon star). All three classes are re-exported from the top-level ``phylozoo``
  namespace.
* New :func:`~phylozoo.core.network.dnetwork.derivations.displayed_triplets` derivation
  on :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`, returning a
  :class:`~phylozoo.core.triplet.tprofileset.TripletProfileSet` aggregated over the
  displayed trees.
* Manual chapter ``Triplets`` (overview + per-class pages) and an ``api/core/triplets``
  API page.
* Manual page ``Parallel Execution`` and ``api/utils/parallel`` API page documenting
  :mod:`phylozoo.utils.parallel` as the standard interface for parallel execution in
  PhyloZoo. No functions currently expose a ``parallel`` parameter; the module is the
  intended way to introduce parallelization in future implementations.
* New ``CI`` GitHub Actions workflow (``.github/workflows/ci.yml``) running ``pytest``
  (Python 3.10 and 3.11), ``ruff``, ``black --check``, and ``mypy`` on every push to
  ``master`` and every pull request. See the :doc:`Testing Guide <testing>` for details.
  A CI badge has been added to the project README.

Changed
^^^^^^^

* :func:`~phylozoo.viz.plot` (and the per-type plotters
  :func:`~phylozoo.viz.dnetwork.plot_dnetwork`,
  :func:`~phylozoo.viz.sdnetwork.plot_sdnetwork`,
  :func:`~phylozoo.viz.d_multigraph.plot_dmgraph`,
  :func:`~phylozoo.viz.m_multigraph.plot_mmgraph`) now have ``show=True`` as the
  default. Pass ``show=False`` to suppress display, for example when saving the figure
  or composing subplots.
* ``ParallelBackend`` class docstring restructured to plain prose (the enum members are
  now described inline rather than in a NumPy-style ``Attributes`` section) to avoid
  duplicate Sphinx object descriptions.
* Black, Ruff, and Mypy fixes have been made thoughout the codebase.

Removed
^^^^^^^

* ``phylozoo.core.distance.operations`` (and the entire ``operations.py`` file). The
  Traveling Salesman Problem solvers (``optimal_tsp_tour``, ``approximate_tsp_tour``)
  have been removed, together with the corresponding manual section and API entry.
  Their implementation is now part of the dependent package ``physquirrel``.
* ``phylozoo.core.quartet.qdistance.quartet_distance_with_partition`` and its
  documentation in the Quartet Profile Sets manual page.

Version 0.1.2
~~~~~~~~~~~~~

* Fixed newick parsing branch length order: it must be after the hybrid marker, not before.

Version 0.1.1
~~~~~~~~~~~~~

* Fixed docs build: added matplotlib to docs dependencies

Version 0.1.0
~~~~~~~~~~~~~

* Initial public release
