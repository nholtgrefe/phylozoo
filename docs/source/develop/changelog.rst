Changelog
=========

Version History
---------------

0.2
~~~

0.2.2
^^^^^

Added
"""""

* ``make_lsa: bool = False`` parameter on
  :func:`~phylozoo.core.network.dnetwork.derivations.displayed_trees`: when ``True``,
  each displayed tree is converted to its LSA network before being returned.
* :attr:`~phylozoo.viz.sdnetwork.style.SDNetStyle.label_rotation` ``(float | None)``
  attribute on :class:`~phylozoo.viz.sdnetwork.style.SDNetStyle`.
  ``None`` (default) auto-aligns each taxon label with the direction of the leaf's
  connecting edge; a fixed float (e.g. ``-7``) applies a uniform tilt to all labels,
  which reduces overlap in dense force-directed layouts.
* ``anchor`` parameter on the internal ``draw_label`` rendering helper
  (:mod:`phylozoo.viz._render`): taxon labels are now automatically placed on the side
  of the node facing away from the connecting edge, preventing labels from crossing edges
  in ``neato`` and related layouts.
* New tutorials:

  * :doc:`../tutorials/visualization_sdnetwork` — step-by-step guide to plotting a
    semi-directed phylogenetic network, using the *Xiphophorus* swordfish example from
    :cite:`Holtgrefe2025Squirrel`.
  * Tutorial on displayed-tree indistinguishable networks: introduces the concept and
    provides worked example code.

Changed
"""""""

* Default layout for :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`
  changed from ``'twopi'`` to ``'neato'`` (Graphviz spring-embedder). All affected call
  sites and documentation updated.
* :class:`~phylozoo.viz.sdnetwork.style.SDNetStyle` default values updated to better
  suit the ``neato`` layout: ``node_color='white'``, ``node_size=80.0``,
  ``leaf_color='#0a0a0a'``, ``leaf_size=100.0``, ``hybrid_color='#fcc0bc'``,
  ``label_offset=0.01`` (previously lighter colours and larger sizes).
* Semi-directed network visualization tutorial and styling manual updated to reflect
  the new style defaults, ``neato`` layout, and ``label_rotation`` parameter.
* Quickstart tutorial revised: network node accessor methods, hybrid network example,
  and plotting instructions updated.
* Citation format in README and documentation updated to a full bibliographic reference
  for :cite:`Holtgrefe2025Squirrel`.

Performance
"""""""""""

* :func:`~phylozoo.core.network.sdnetwork.derivations.k_taxon_subnetworks`: all
  :math:`\binom{n}{2}` pairwise up-down vertex sets are now pre-computed once and
  cached before the combination loop.  Each :func:`~phylozoo.core.network.sdnetwork.derivations.subnetwork`
  call performs a dictionary lookup instead of re-invoking
  :func:`~phylozoo.core.primitives.m_multigraph.features.updown_path_vertices`.
  For :math:`k=5` on :math:`n=13` taxa this reduces ``updown_path_vertices`` calls from
  :math:`\binom{13}{5} \cdot \binom{5}{2} = 12{,}870` to :math:`\binom{13}{2} = 78`, a
  :math:`{\approx}165\times` reduction.
* :func:`~phylozoo.core.network.sdnetwork.derivations.subnetwork`: removed a redundant
  ``MixedMultiGraph.copy()`` call.
  :func:`~phylozoo.core.primitives.m_multigraph.transformations.subgraph` already
  returns a fresh independent object; the subsequent ``.copy()`` was an unnecessary
  second allocation.
* :func:`~phylozoo.core.network.sdnetwork.transformations.identify_parallel_edges` /
  :func:`~phylozoo.core.network.sdnetwork.derivations.subnetwork`: extracted
  ``_identify_parallel_edges_inplace(graph, exclude_nodes)`` as an internal helper that
  modifies a :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`
  in-place.  ``subnetwork`` now calls this helper directly on the working graph before
  the single :func:`~phylozoo.core.network.sdnetwork.conversions.sdnetwork_from_graph`
  call, eliminating one ``MixedMultiGraph.copy()`` and one
  :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`
  construction per call when ``identify_parallel_edges=True``.  The public
  ``identify_parallel_edges`` function is an unchanged thin wrapper.
* **Combined effect**: on a level-3 network with :math:`k=11` backbone leaves
  (14 taxa, 2,002 quinnets) quinnet extraction time fell from ≈22 s to ≈4 s
  (≈5× end-to-end speedup).

Fixed
"""""

* License badge in ``README.md`` corrected to display MIT.

Internal
""""""""

* Distance matrix helper functions in
  :mod:`phylozoo.core.distance.classifications` and
  :mod:`phylozoo.core.distance.decomposition` refactored for readability;
  no API change.
* Edge routing for the ``pz-radial`` and spring-embedder (``neato``/``fdp``) layouts
  improved in :mod:`phylozoo.viz.sdnetwork.layout`.

0.2.1
^^^^^

Added
"""""

* Split decomposition (:cite:`Bandelt1992`) for distance matrices:

  * :func:`~phylozoo.core.distance.decomposition.isolation_index` — computes the isolation
    index of a bipartition with respect to a distance matrix.
  * :func:`~phylozoo.core.distance.decomposition.split_decomposition` — canonical
    decomposition ``d = d^0 + Σ α_S δ_S``, returning a
    :class:`~phylozoo.core.split.weighted_splitsystem.WeightedSplitSystem` of all d-splits
    and the split-prime residual as a :class:`~phylozoo.core.distance.base.DistanceMatrix`.
  * Both functions are re-exported from ``phylozoo.core.distance``.

* New classification functions in :mod:`phylozoo.core.distance.classifications`:

  * :func:`~phylozoo.core.distance.classifications.is_tree_metric` — four-point condition
    check (Numba-accelerated, :math:`O(n^4)`).
  * :func:`~phylozoo.core.distance.classifications.is_totally_decomposable` — checks
    whether the split-prime residual is zero.

* :func:`~phylozoo.core.split.algorithms.distances_from_splitsystem` ↔
  :func:`~phylozoo.core.distance.decomposition.split_decomposition` round-trip tests and
  :func:`~phylozoo.core.split.algorithms.tree_from_splitsystem` ↔
  :func:`~phylozoo.core.network.sdnetwork.derivations.induced_splits` round-trip tests.
* Documentation: new *Split Decomposition* section in the distance manual; updated
  *Algorithms* section in the split-system manual; new ``api/core/distance`` entry for
  the decomposition module. BibTeX entry ``Bandelt1992`` added to ``bibliography.bib``.

Changed
"""""""

* ``_check_kalmanson_conditions`` moved from a nested closure inside
  :func:`~phylozoo.core.distance.classifications.is_kalmanson` to a module-level
  ``@njit`` function, consistent with the other Numba helpers in the module.

0.2.0
^^^^^

Added
"""""

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
  (Python 3.10 and 3.11), ``ruff``, and ``black --check`` on every push to ``master``
  and every pull request. See the :doc:`Testing Guide <testing>` for details.
  A CI badge has been added to the project README.

Changed
"""""""

* ``ParallelBackend`` class docstring restructured to plain prose (the enum members are
  now described inline rather than in a NumPy-style ``Attributes`` section) to avoid
  duplicate Sphinx object descriptions.
* Black, Ruff, and Mypy fixes have been made throughout the codebase.

Removed
"""""""

* ``phylozoo.core.distance.operations`` (and the entire ``operations.py`` file). The
  Traveling Salesman Problem solvers (``optimal_tsp_tour``, ``approximate_tsp_tour``)
  have been removed, together with the corresponding manual section and API entry.
  Their implementation is now part of the dependent package ``physquirrel``.
* ``phylozoo.core.quartet.qdistance.quartet_distance_with_partition`` and its
  documentation in the Quartet Profile Sets manual page.

0.1
~~~

0.1.2
^^^^^

* Fixed newick parsing branch length order: it must be after the hybrid marker, not before.

0.1.1
^^^^^

* Fixed docs build: added matplotlib to docs dependencies.

0.1.0
^^^^^

* Initial public release.
