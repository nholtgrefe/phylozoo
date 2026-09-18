Plotting
========

The :mod:`phylozoo.viz` module provides plotting for phylogenetic networks and graphs.
Plots are rendered using `matplotlib <https://matplotlib.org/stable/>`_ and return
:class:`matplotlib.axes.Axes` objects.

The :func:`~phylozoo.viz.plot` function accepts any of four object types and dispatches by type to the
appropriate plotter:

* :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` → :func:`~phylozoo.viz.dnetwork.plot_dnetwork`
* :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` → :func:`~phylozoo.viz.sdnetwork.plot_sdnetwork`
* :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` → :func:`~phylozoo.viz.d_multigraph.plot_dmgraph`
* :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` → :func:`~phylozoo.viz.m_multigraph.plot_mmgraph`

How to Plot
-----------

The :func:`~phylozoo.viz.plot` function has the following parameters:

* **obj** — The object to plot (one of the four supported types).
* **layout** (str, default='auto') — Layout algorithm. Use ``'auto'`` for the default per type, or specify a name (e.g. ``'pz-cladogram'``, ``'spring'``, ``'circular'``). See the layout section below for more details.
* **style** — Style object. See the :doc:`Styling <styling>` documentation for more details. Use ``None`` for the default style for the object type.
* **ax** — Optional matplotlib axes. Plot on existing axes (e.g. for subplots).
* **show** (bool, default=False) — If ``True``, display the plot. If ``False``, return the axes for saving or customization.
* **\*\*kwargs** — Layout-specific parameters (e.g. ``layer_gap=2.0`` for ``'pz-cladogram'``).

For quickly plotting a network, and you don't need to customize the layout or style, you can use the following simple code to plot the network:

.. code-block:: python

   from phylozoo.viz import plot
   from phylozoo import DirectedPhyNetwork

   dnet = DirectedPhyNetwork.load("network.enewick")
   plot(dnet)

Layouts
-------

.. _viz-layout:

A layout algorithm determines where nodes are placed. Pass its name as the ``layout`` argument;
layout-specific parameters go in ``**kwargs``:

.. code-block:: python

   plot(dnet)                                                # default: pz-cladogram
   plot(dnet, layout='pz-cladogram', rectangular=True)
   plot(dnet, layout='pz-layered', direction='LR')
   plot(sdnet)                                               # default: pz-unrooted
   plot(sdnet, layout='pz-unrooted', daylight=8, refine=0)
   plot(sdnet, layout='pz-radial', root_location=some_node)

PhyloZoo layouts
^^^^^^^^^^^^^^^^

Four layouts designed for phylogenetic networks, all deterministic and needing nothing beyond
``phylozoo[viz]``: two rooted, layered ones for directed networks and two unrooted-style ones available
for both classes. These are the recommended layouts and the defaults.

* **pz-cladogram** (:class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` only; default) — Layered drawing
  with a tree backbone. Every node hangs below its lowest parent; the other parent edges of hybrid nodes are drawn as
  reticulate edges. Siblings are ordered by a hybrid-aware barycenter heuristic followed by a local search on the
  exact number of edge crossings, so sibling subtrees stay contiguous and reticulate edges stay short. Arrowheads are
  drawn on hybrid edges only (see :attr:`~phylozoo.viz.dnetwork.style.DNetStyle.arrows`).
  Parameters:

  - ``layer_gap`` (float, default 1.0) — Spacing between layers;
  - ``leaf_gap`` (float, default 1.0) — Spacing between consecutive leaves;
  - ``trials`` (int, default 5) — Number of ordering attempts (the first from the node order, the rest from random orders);
  - ``seed`` (int or None, default 0) — Random seed for the restarts;
  - ``direction`` (str, default ``'TD'``) — ``'TD'`` (root at the top) or ``'LR'`` (root on the left);
  - ``align_leaves`` (bool, default True) — Put all leaves on the bottom layer;
  - ``rectangular`` (bool, default False) — Draw edges as orthogonal elbows (backbone edges along the parent's layer
    and then into the child, reticulate edges entering the hybrid sideways along its layer);
  - ``horizontal_reticulations`` (bool or None, default None) — Lower the parent of every reticulate edge onto the
    hybrid's layer where the network allows it, so that edge is a single horizontal segment. ``None`` means on for
    ``rectangular`` drawings and off otherwise;
  - ``x_scale`` / ``y_scale`` (float, default 1.0) — Scaling factors for the coordinates.

  Recommended: ``plot(net)`` for a cladogram-style picture, ``plot(net, rectangular=True)`` for the classic
  rectangular network drawing, ``plot(net, direction='LR')`` with horizontal leaf labels for long taxon names.
* **pz-layered** (:class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` only) — Self-contained layered
  (Sugiyama) layout in the spirit of Graphviz ``dot``: longest-path layers, edges spanning several layers routed
  through bend points, layer-by-layer barycenter sweeps with adjacent swaps to reduce crossings, and a coordinate
  pass that straightens edges. Nodes of a layer are ordered freely (no tree backbone), which gives fewer crossings on
  heavily reticulate networks than ``pz-cladogram`` at the cost of sibling subtrees not staying together; leaves keep their
  own depth by default.
  Parameters:

  - ``layer_gap`` (float, default 1.0) — Spacing between layers;
  - ``node_gap`` (float, default 1.0) — Minimum spacing between nodes of a layer;
  - ``sweeps`` (int, default 8) — Number of down-and-up ordering sweeps;
  - ``direction`` (str, default ``'TD'``) — ``'TD'`` or ``'LR'``;
  - ``align_leaves`` (bool, default False) — Put all leaves on the bottom layer;
  - ``rectangular`` (bool, default False) — Draw every edge segment as an orthogonal elbow; the last segment into a
    hybrid node arrives sideways;
  - ``x_scale`` / ``y_scale`` (float, default 1.0) — Scaling factors for the coordinates.

  Recommended for networks with many reticulations where ``pz-cladogram`` shows long crossing reticulate edges.
* **pz-radial** (both classes) — Circular cladogram, as in the radial network view of Dendroscope: root at the
  centre, leaves evenly spaced on the outer circle, a node's radius given by its layer and its angle by the mean of
  its children. The tree backbone and the leaf order come from the ``pz-cladogram`` computation, so reticulate edges
  (straight chords) are kept short. A semi-directed network is rooted first with
  :func:`~phylozoo.core.network.sdnetwork.derivations.to_d_network`.
  Parameters:

  - ``radius`` (float, default 1.0) — Radius of the leaf circle;
  - ``start_angle`` (float, default 0.0) — Angle of the first leaf, in radians;
  - ``angle_direction`` (str, default ``'clockwise'``) — ``'clockwise'`` or ``'counterclockwise'``;
  - ``root_location`` (semi-directed only; node or edge, default None) — Where to root the network;
  - ``trials`` / ``seed`` — Ordering options of ``pz-cladogram``.
* **pz-unrooted** (both classes; default for
  :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`) — Unrooted tree-of-blobs
  layout. The network is decomposed into blobs; the tree obtained by contracting them is drawn with the equal-angle
  algorithm (rooted at its centroid), and each blob is drawn as a regular polygon whose outer cycle comes from a
  planar embedding, with interior nodes at the barycentre of their neighbours. Pendant subtrees leave a blob in the
  cyclic order of their attachment nodes, so cut edges never cross blob edges; for trees this is the classic
  equal-angle unrooted tree drawing. Two optional finishing steps spread the drawing: ``daylight`` rounds rotate
  the rigid subtrees around every tree node and around every blob until the free angles between them are equal
  (every edge keeps its exact length), and ``refine`` iterations of stress majorization (the ``neato`` algorithm,
  started from this crossing-free drawing) relax the whole picture; the stress step only ever accepts a state with no
  more edge crossings than the drawing it started from. On a directed network the root is an ordinary
  node and every edge carries an arrowhead.
  Parameters:

  - ``edge_length`` (float, default 1.0) — Length of cut edges (between blobs and to leaves);
  - ``node_spacing`` (float, default 0.7) — Target distance between consecutive nodes on a blob's circle;
  - ``inner_scale`` (float, default 0.6) — Scale of subtrees that hang off interior blob nodes and are drawn inside the blob;
  - ``start_angle`` (float, default 0.0) — Rotation of the whole drawing, in radians;
  - ``daylight`` (int, default 0) — Equal-daylight rounds, applied first;
  - ``refine`` (int, default 50) — Stress-majorization iterations, applied last (``0`` disables it).

  Recommended settings: the default (stress only) for a compact ``neato``-like picture with round blobs;
  ``daylight=8, refine=0`` for the SplitsTree-style unrooted drawing with exact edge lengths;
  ``daylight=4`` together with the default ``refine`` when the stress step folds subtrees on top of each other;
  ``node_spacing=1.0`` for uniform edge lengths (plain ``neato`` look); ``refine=0, daylight=0`` for the raw
  geometry with blobs as exact polygons.

Other layouts
^^^^^^^^^^^^^

The generic graph layouts of NetworkX (``spring``, ``kamada_kawai``, ``circular``, ``shell``,
``spectral``, ``spiral``, ``planar``, ``random``, ``bipartite``) can be passed as ``layout`` for all four
object types; their keyword arguments (``k``, ``iterations``, ``seed``, ...) go in ``**kwargs``. If
`PyGraphviz <https://pygraphviz.github.io/>`_ is installed (see :doc:`Installation <../installation>`),
the Graphviz programs ``dot``, ``neato``, ``twopi``, ``fdp``, ``sfdp`` and ``circo`` are available the same
way, with Graphviz attributes passed through ``args``. Several PhyloZoo layouts build on the same ideas
(``pz-layered`` on ``dot``, ``pz-radial`` on ``twopi``, the refinement of ``pz-unrooted`` on ``neato``) but
are tuned to phylogenetic networks, so for networks the PhyloZoo layouts are the better choice. The generic
layouts are mainly there for plotting :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`
and :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`, for which they are the only
layouts, with ``spring`` as the default.

Text Drawings
-------------

:func:`~phylozoo.viz.to_pretty_print` (also available as the ``to_pretty_print()`` and ``pretty_print()`` methods of both
network classes) draws a network as text, root on the left and leaves on the right, without matplotlib:

.. code-block:: python

   >>> net.pretty_print()
   ┌───────────────●───────● C
   │               v
   ○       ┊┄┄┄┄┄┄>◆───────● A
   └───────●
           └───────────────● B

Tree edges are box-drawing lines, hybrid edges dotted lines ending in an arrowhead at the hybrid node ``◆``;
``○`` is the root, ``●`` any other node; ``┼`` marks a hybrid edge crossing a tree edge. The placement is that of
``pz-cladogram`` with ``rectangular=True``, so hybrid edges are horizontal wherever the network allows it and sibling
subtrees stay together. Parallel edges are drawn as a loop (the second copy detours one row). Semi-directed
networks are rooted first (``root_location`` chooses where).

Large networks are compacted automatically (narrower columns down to 3 characters until the drawing fits
``max_width``, default 120, and no blank rows between leaves above 25 leaves); networks with more than
``max_leaves`` (default 80) leaves are refused, since one row per leaf stops being readable. Use ``col_width`` and
``rows_per_leaf`` to override the automatic choices.

Saving a Figure
---------------

The :func:`~phylozoo.viz.plot` function returns a matplotlib axes object. Set ``show=False`` and call
the matplotlib ``savefig`` method on the figure:

.. code-block:: python

   ax = plot(network, show=False)
   ax.figure.savefig("network.png", dpi=300, bbox_inches="tight")

Here, ``dpi`` controls the resolution of the saved image and ``bbox_inches='tight'`` crops the figure to the smallest bounding box that contains all plot elements, removing excess whitespace and ensuring labels are not cut off.


See Also
--------

- :doc:`Styling <styling>` — Colors, sizes, and appearance
- :doc:`Tutorial: Plotting a Semi-Directed Network <../../tutorials/visualization_sdnetwork>` — End-to-end example with layout and style customisation
