Styling
=======

The :mod:`phylozoo.viz` module provides style classes for controlling the visual
appearance of network plots: colors, sizes, and labels.

Using Styles
------------

Each of the four plottable classes has its own style class:

* :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork` uses :class:`~phylozoo.viz.dnetwork.style.DNetStyle`
* :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork` uses :class:`~phylozoo.viz.sdnetwork.style.SDNetStyle`
* :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph` uses :class:`~phylozoo.viz.d_multigraph.style.DMGraphStyle`
* :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph` uses :class:`~phylozoo.viz.m_multigraph.style.MGraphStyle`

Styles inherit from a base hierarchy: :class:`~phylozoo.viz.d_multigraph.style.BaseStyle` is the root; :class:`~phylozoo.viz.d_multigraph.style.DMGraphStyle` and :class:`~phylozoo.viz.m_multigraph.style.MGraphStyle` extend it for
the two graph types; ``DNetStyle`` extends ``DMGraphStyle`` with leaf and hybrid options for directed networks;
``SDNetStyle`` extends ``MGraphStyle`` with the same options for semi-directed networks.

Import the style from the corresponding viz submodule:

.. code-block:: python

   from phylozoo.viz.dnetwork import DNetStyle
   from phylozoo.viz.sdnetwork import SDNetStyle
   from phylozoo.viz.d_multigraph import DMGraphStyle
   from phylozoo.viz.m_multigraph import MGraphStyle

Each submodule provides a ``default_style()`` function that returns a style instance with default values.
If you omit the ``style`` argument when calling :func:`~phylozoo.viz.plot`, the default style for that object type
is used automatically.

Styling Parameters
------------------

BaseStyle
^^^^^^^^^

The :class:`~phylozoo.viz.d_multigraph.style.BaseStyle` class is the base class for all styles. It defines the common attributes for all styles.

* **node_color** (str, default='lightblue'): Color for nodes. Any valid matplotlib color string is supported.
* **node_size** (float, default=500.0): Size of nodes.
* **edge_color** (str, default='gray'): Color for edges. Any valid matplotlib color string is supported.
* **edge_width** (float, default=2.0): Width of edges.
* **with_labels** (bool, default=True): Whether to show node labels.
* **label_offset** (float, default=0.12): Offset for labels from nodes.
* **label_font_size** (float, default=10.0): Font size for labels.
* **label_color** (str, default='black'): Color for labels. Any valid matplotlib color string is supported.

DMGraphStyle and MGraphStyle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :class:`~phylozoo.viz.d_multigraph.style.DMGraphStyle` class extends :class:`~phylozoo.viz.d_multigraph.style.BaseStyle` for :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`.
The :class:`~phylozoo.viz.m_multigraph.style.MGraphStyle` class extends :class:`~phylozoo.viz.d_multigraph.style.BaseStyle` for :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`.
The two classes hold the exact same attributes as :class:`~phylozoo.viz.d_multigraph.style.BaseStyle`, but are defined with future extensibility in mind.

DNetStyle
^^^^^^^^^

Extends DMGraphStyle for :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`.
Adds phylogenetic network-specific options:

* **leaf_color** (str, default='lightblue'): Color for leaf nodes. Any valid matplotlib color string is supported.
* **hybrid_color** (str, default='lightblue'): Color for hybrid nodes. Any valid matplotlib color string is supported.
* **leaf_size** (float | None, default=None): Size of leaf nodes; if None, uses node_size
* **hybrid_edge_color** (str, default='red'): Color for hybrid edges. Any valid matplotlib color string is supported.

SDNetStyle
^^^^^^^^^^

Extends MGraphStyle for :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`.
Overrides several base defaults to suit semi-directed network plots and adds
phylogenetic options:

* **node_color** (str, default='white'): Color for internal nodes.
* **node_size** (float, default=80.0): Size of internal nodes.
* **label_offset** (float, default=0.01): Label distance from node centre in layout coordinates.
* **leaf_color** (str, default='#0a0a0a'): Color for leaf nodes.
* **leaf_size** (float, default=100.0): Size of leaf nodes.
* **hybrid_color** (str, default='#fcc0bc'): Color for hybrid nodes.
* **hybrid_edge_color** (str, default='red'): Color for hybrid edges.
* **label_rotation** (float | None, default=None): Rotation of taxon labels in degrees. ``None`` aligns each label with its pendant edge (reading outwards); a fixed value such as ``0`` keeps all labels at that angle and places them beside the node on its outward side.
  ``None`` (the default) auto-aligns each label with the direction of its connecting edge,
  which reduces overlap in force-directed layouts like ``neato``.
  A fixed float (e.g. ``-7``) applies a uniform tilt to all labels.

Example
-------

The following example loads a network, creates a custom style with several parameters,
and plots the result:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.viz import plot
   from phylozoo.viz.dnetwork import DNetStyle

   network = DirectedPhyNetwork.load("network.enewick")

   style = DNetStyle(
       node_color="steelblue",
       leaf_color="lightgreen",
       hybrid_color="salmon",
       hybrid_edge_color="red",
       node_size=150,
       leaf_size=200,
       edge_color="darkgray",
       edge_width=2.5,
       with_labels=True,
       label_font_size=12,
   )

   plot(network, style=style)

See Also
--------

- :doc:`Plotting <plotting>` — How to plot with a style
- :doc:`viz API reference <../../api/viz/index>` — Full style class documentation
