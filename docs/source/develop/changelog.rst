Changelog
=========

Version History
---------------

0.3
~~~

0.3.0
^^^^^

A performance release: the fundamental operations were profiled on networks of up
to :math:`10^5` leaves and the algorithmic bottlenecks removed, without changing
results. Every docstring example is now executed as a test.

Added
"""""

* ``copy`` keyword on :func:`~phylozoo.core.network.dnetwork.conversions.dnetwork_from_graph`
  and :func:`~phylozoo.core.network.sdnetwork.conversions.sdnetwork_from_graph`.
  With ``copy=False`` the network adopts the given graph instead of rebuilding it
  (the caller must not use the graph afterwards); label validation, leaf
  auto-labelling and :meth:`validate` run exactly as before. All internal
  transformations use it, which removes one full graph construction from every
  conversion, subnetwork and displayed tree.
* :meth:`has_node` and :meth:`all_degrees` on
  :class:`~phylozoo.core.primitives.d_multigraph.DirectedMultiGraph` and
  :class:`~phylozoo.core.primitives.m_multigraph.MixedMultiGraph` (``all_degrees``
  returns every node's in-/out-/undirected degree in one pass), and
  :func:`~phylozoo.core.primitives.m_multigraph.features.has_updown_path` (linear-time
  existence test for an up-down path, without enumerating paths).
* Docstring examples run as doctests in the normal ``pytest`` run
  (``--doctest-modules`` and ``src/phylozoo`` in ``testpaths``); see the
  :doc:`Testing guide <testing>`. About 120 stale examples were repaired, and
  examples that wrote files now write into a temporary directory.
* README quickstart, and a manual section on
  :ref:`validation inside library functions <validation-inside-library-functions>`.

Changed
"""""""

* **Derived networks are no longer re-validated.** Functions that build a new
  network from an already valid one — ``subnetwork``, ``to_sd_network``,
  ``to_d_network`` (with an automatically chosen root), ``tree_of_blobs``,
  ``displayed_trees``, ``to_lsa_network``, ``suppress_2_blobs``,
  ``identify_parallel_edges``, ``binary_resolution`` and the generator enumeration —
  construct their result inside :func:`~phylozoo.utils.validation.no_validation`,
  since it is valid by construction. Constructors, parsers and functions whose result
  depends on caller input (e.g. ``to_d_network`` with an explicit ``root_location``,
  ``root_at_outgroup``) still validate. Roughly halves the default cost of
  ``to_sd_network``, ``displayed_trees`` and ``suppress_2_blobs``.
* :func:`~phylozoo.core.network.dnetwork.derivations.distances` and
  :func:`~phylozoo.core.network.dnetwork.derivations.displayed_splits` (both
  representations) no longer enumerate all :math:`\prod_b 2^{r_b}` switchings. A
  path's stretch inside a blob depends only on that blob's hybrid choices, so only a
  reference switching plus each blob's local switchings are evaluated
  (:math:`1 + \sum_b 2^{r_b}`) and aggregated per blob — for all three distance modes,
  since a sum of independently chosen terms is minimised or maximised term by term.
  The cost is exponential in the *level*, not in the number of reticulations: 20
  reticulations in 20 blobs take 1.6 s instead of :math:`2^{20}` switchings. Both
  functions share one switching machinery that walks the network's cached graph with
  removed edges instead of copying a switching graph per switching (single-blob
  ``distances`` 2.7x faster; the semi-directed matrix computation, which used one
  shortest-path call per ordered leaf pair, is ~20x faster).
* :func:`~phylozoo.core.distance.classifications.is_tree_metric` is
  :math:`O(n^2)` instead of :math:`O(n^4)`: the tree an additive matrix must come from
  is reconstructed by inserting taxa one by one, and the matrix is a tree metric
  exactly when that tree reproduces it (n=500: 4.2 s to 0.26 s).
  :func:`~phylozoo.core.split.algorithms.tree_from_splitsystem` inserts compatible
  splits as clusters, smallest first, instead of searching cut-vertices per split
  (n=120: 2.1 s to 0.025 s).
  :func:`~phylozoo.core.sequence.distances.hamming_distances` counts matches and valid
  sites with (exact, blocked) matrix products instead of a pair loop (6–31x);
  :func:`~phylozoo.core.split.classifications.is_pairwise_compatible` uses bitmasks.
* Displayed trees / quartets / triplets: no per-tree graph copies, one shared rooting
  for semi-directed networks, degree-1 pruning as a worklist, and quartets/triplets
  read topology off the tree graphs instead of building network objects.
* Semi-directed validation no longer round-trips through a directed network;
  ``lsa_node`` uses dominators; ``subnetwork`` collects ancestors in one traversal;
  ``is_ultrametric``, ``is_normal``, ``source_components``, ``orient_away_from_vertex``
  and ``quartet_distance`` avoid repeated scans.
* Primitives: the combined undirected view and the ``nodes``/``edges`` views are built
  lazily (constructing a network no longer maintains a second graph), ``subgraph``
  walks only the selected nodes' adjacency (extracting all blobs of a 2 000-leaf
  network: 1.37 s to 0.03 s), and ``Split`` stores its sides as the partition's own
  frozensets and can share one taxon set across a split system (2.9x less memory for
  large split systems; ``split.set1``/``set2`` are now ``frozenset`` — equality with
  plain sets is unchanged, and so is the ``repr``).
* Generators: R1 is applied to unordered side pairs (it is symmetric), the reachability
  matrix uses set unions, and ``semidirect_generators`` no longer re-validates its
  by-construction-valid results (level-4 semi-direction 8.2 s to 1.2 s). Generator
  counts are unchanged.

Fixed
"""""

* ``distances()`` raised ``PhyloZooValueError: Multiple parallel edges exist`` on any
  network with parallel edges (branch lengths were looked up without the edge key).
* ``displayed_quartets`` on semi-directed networks crashed because a switching kept
  the surviving hybrid edges directed; switchings are now undirected before use.
* The DOT writers dropped non-zero edge keys on non-parallel edges, so such graphs did
  not round-trip (the cause of an intermittently failing generator I/O test).
* Nested :func:`~phylozoo.utils.validation.no_validation` blocks now stack as
  documented: an inner block adds to the outer block's suppression instead of
  replacing it when it names other methods.
* Documentation: duplicate Squirrel citation removed.

0.2
~~~

0.2.6
^^^^^

Added
"""""

* Standard Graphviz ``dot`` format for
  :class:`~phylozoo.core.primitives.m_multigraph.MixedMultiGraph`, alongside the
  existing ``phylozoo-dot`` (which remains the default). Undirected edges are
  written as ``u -> v [dir=none]`` inside a ``digraph``, so — unlike
  ``phylozoo-dot``, which mixes ``->`` and ``--`` in one ``graph`` block and is not
  valid DOT — the output opens in any Graphviz tool and round-trips losslessly,
  including parallel directed and undirected edges. The DOT scaffolding (escaping,
  attribute formatting/parsing, node-id coercion and document parsing) is now
  shared between the directed- and mixed-multigraph handlers in a new
  :mod:`phylozoo.utils.io.format_utils.dot` module.
* :class:`~phylozoo.core.network.sdnetwork.generator.SemiDirectedGenerator` and
  :class:`~phylozoo.core.network.dnetwork.generator.DirectedGenerator` now inherit
  :class:`~phylozoo.utils.io.IOMixin`, gaining ``save``/``load``/``to_string``/
  ``from_string``/``convert``. They inherit the formats of their underlying graph
  — ``phylozoo-dot`` (default) and ``dot`` for the semi-directed generator,
  ``dot`` (default) and ``edgelist`` for the directed one — so a generator can be
  serialised and reconstructed directly, without routing through
  ``generator.graph``. Reading rebuilds (and validates) the generator from the
  parsed graph.
* Generator construction is now exposed as reusable single steps:
  :func:`~phylozoo.core.network.dnetwork.generator.construction.gambette_step`
  (apply the Gambette R1/R2 rules + isomorphism deletion to a collection of
  directed generators) and
  :func:`~phylozoo.core.network.sdnetwork.generator.construction.semidirect_generators`
  (semi-direct a collection of directed generators + isomorphism deletion). The
  ``all_level_k_generators`` functions are now thin loops over these, and a
  caller can advance / semi-direct their own (e.g. saved) generator set without
  rebuilding from level 0. Both accept mixed-level inputs (``gambette_step`` also
  accepts level-0 generators) and assume their inputs are valid.

Changed
"""""""

* Generator enumeration now deduplicates with a Weisfeiler-Lehman graph hash added
  to the cheap-invariant key (new helpers ``_get_graph_wl_hash`` on the directed
  and mixed multigraph isomorphism modules; multiplicity- and direction-aware via a
  simple-graph encoding). This splits isomorphism candidates into far finer groups,
  so the exact VF2 check runs only within tiny groups — the level-4 directed
  enumeration drops from ~276 s to ~6 s (~44x), and the full level-4 semi-directed
  enumeration from ~285 s to ~9 s. Output is identical (1993 directed / 307
  semi-directed level-4 generators); exactness is unchanged, since the hash is only
  a grouping key and VF2 remains the final check.

Fixed
"""""

* ``phylozoo-dot`` round-trip for
  :class:`~phylozoo.core.primitives.m_multigraph.MixedMultiGraph`: the writer
  injected a ``label=<node id>`` attribute on every label-less node, and the reader
  parsed that unquoted value back as an ``int``. A round-tripped graph therefore
  carried a spurious **non-string** ``label`` on every node, which broke downstream
  label validation — e.g.
  :func:`~phylozoo.core.network.sdnetwork.generator.attachment.attach_leaves_to_generator`
  raised ``PhyloZooTypeError: ... non-string label``. The writer no longer emits the
  redundant label, and the reader keeps any ``label`` value as a string. Generator
  graphs (and any ``MixedMultiGraph``) now round-trip losslessly, including parallel
  directed and undirected edges, and a loaded generator is directly usable to build
  networks — without first rebuilding a clean graph from its edge lists.

0.2.5
^^^^^

Fixed
"""""

* :attr:`~phylozoo.core.network.sdnetwork.generator.SemiDirectedGenerator.hybrid_sides`:
  a reticulation whose child slot was already occupied by an undirected (backbone)
  edge was still reported as a hybrid *side*. Attaching the required leaf to it then
  produced an invalid hybrid (total degree = in-degree + 2), raising a degree error
  in :func:`~phylozoo.core.network.sdnetwork.generator.attachment.attach_leaves_to_generator`.
  ``hybrid_sides`` now only includes hybrid nodes whose child slot is free (no
  directed out-edge **and** no incident undirected edge); such occupied
  reticulations remain hybrid *nodes* (the generator's level is unchanged) and
  their descendants are reached through the incident undirected edge side. The
  level-1 bidirected self-loop node is still a hybrid side.

0.2.4
^^^^^

Fixed
"""""

* :func:`~phylozoo.core.network.sdnetwork.derivations.subnetwork`: the updown-path
  cache lookup used an ``or``-chain (``cache.get(key) or ...``), which treated an
  empty cached path set (falsy) as a cache miss and silently recomputed the path.
  Fixed by using an explicit ``_MISS_SENTINEL = object()`` sentinel so empty sets are
  handled correctly.

0.2.3
^^^^^

* Fix PyPi deployment

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
