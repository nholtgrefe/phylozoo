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
* **layout** (str, default='auto') — Layout algorithm. Use ``'auto'`` for the default per type, or specify a name (e.g. ``'pz-dag'``, ``'spring'``, ``'circular'``). See the layout section below for more details.
* **style** — Style object. See the :doc:`Styling <styling>` documentation for more details. Use ``None`` for the default style for the object type.
* **ax** — Optional matplotlib axes. Plot on existing axes (e.g. for subplots).
* **show** (bool, default=False) — If ``True``, display the plot. If ``False``, return the axes for saving or customization.
* **\*\*kwargs** — Layout-specific parameters (e.g. ``layer_gap=2.0`` for ``'pz-dag'``).

For quickly plotting a network, and you don't need to customize the layout or style, you can use the following simple code to plot the network:

.. code-block:: python

   from phylozoo.viz import plot
   from phylozoo import DirectedPhyNetwork

   dnet = DirectedPhyNetwork.load("network.enewick")
   plot(dnet)

Layouts
-------

.. _viz-layout:

A layout algorithm determines how nodes are positioned on the canvas. PhyloZoo supports
three sources of layouts: **NetworkX** (built-in, no extra dependencies beyond ``phylozoo[viz]``),
**Graphviz** (requires ``pip install phylozoo[graphviz]`` and the Graphviz system library),
and **PhyloZoo** (custom layouts for directed and semi-directed networks).

Pass the layout name as the ``layout`` argument; layout-specific parameters (e.g. ``layer_gap``
for ``pz-dag``) go in ``**kwargs``:

.. code-block:: python

   plot(network, layout='spring', k=2.0, iterations=50, show=True)
   plot(network, layout='pz-dag', layer_gap=2.0, leaf_gap=1.5, show=True)

NetworkX and Graphviz layouts are supported by all four object types. PhyloZoo layouts
are type-specific (see below).

NetworkX
^^^^^^^^

NetworkX layouts are built-in and require no extra installation. Common parameters:
``k`` (optimal distance between nodes), ``iterations``, ``seed``.

* **spring** — Force-directed layout (Fruchterman-Reingold). Well-suited for general networks and graphs.
* **circular** — Nodes arranged on a circle. Good for compact overviews.
* **kamada_kawai** — Force-directed layout. Often produces evenly spaced, readable layouts for networks.
* **planar** — Planar embedding (only for planar graphs).
* **random** — Random node positions. Useful for testing.
* **shell** — Nodes in concentric shells.
* **spectral** — Uses eigenvectors of the graph Laplacian.
* **spiral** — Spiral arrangement of nodes.
* **bipartite** — Two-column layout for bipartite graphs.

For full parameter documentation (e.g. ``k``, ``iterations``, ``seed``, ``scale``, ``center``), see the
`NetworkX layout reference <https://networkx.org/documentation/stable/reference/drawing.html#graph-layout>`_.

Graphviz
^^^^^^^^

Graphviz layouts produce high-quality layouts and scale well to larger graphs. Install with:

.. code-block:: bash

   pip install phylozoo[graphviz]

You must also install the Graphviz system library (e.g. ``apt install graphviz graphviz-dev``
on Debian/Ubuntu, ``brew install graphviz`` on macOS). See the
`PyGraphviz installation guide <https://pygraphviz.github.io/documentation/stable/install.html>`_
for details.

* **dot** — Hierarchical (layered) layout. Well-suited for directed acyclic graphs and phylogenetic trees.
* **twopi** — Radial layout with root at center. Good for trees and semi-directed networks (default for SemiDirectedPhyNetwork).
* **neato** — Spring-model layout.
* **fdp** — Force-directed placement.
* **sfdp** — Scalable force-directed layout for large graphs.
* **circo** — Circular layout.

For layout program parameters and attributes, see the
`PyGraphviz documentation <https://pygraphviz.github.io/documentation/stable/>`_ and the
`Graphviz documentation <https://graphviz.org/documentation/>`_.


PhyloZoo
^^^^^^^^

Custom layouts optimized for phylogenetic networks. Only available for specific types:

* **pz-dag** (:class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` only) — Tree-backbone layout with crossing minimization.
  Parameters:

  - ``layer_gap`` (float, default 1.5) — Spacing between hierarchical layers;
  - ``leaf_gap`` (float, default 1.0) — Spacing between leaves within a layer;
  - ``trials`` (int, default 2000) — Number of random child orderings to try for crossing minimization;
  - ``seed`` (int or None) — Random seed for reproducibility;
  - ``direction`` (str, default ``'TD'``) — Layout direction: ``'TD'`` (top-down) or ``'LR'`` (left-right);
  - ``x_scale`` (float, default 1.5) — Scaling factor for x coordinates;
  - ``y_scale`` (float, default 1.0) — Scaling factor for y coordinates.
* **pz-radial** (:class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` only, trees only) — Radial layout with root at center.
  Parameters:

  - ``radius`` (float, default 1.0) — Maximum radius for leaf nodes;
  - ``start_angle`` (float, default 0.0) — Starting angle in radians for the first leaf;
  - ``angle_direction`` (str, default ``'clockwise'``) — Angle progression: ``'clockwise'`` or ``'counterclockwise'``.

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
