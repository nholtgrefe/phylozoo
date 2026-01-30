Styling
=======

Styling controls the visual appearance of network plots, including colors, sizes, 
and labels. PhyloZoo provides style classes for each network and graph type.

Style Classes
-------------

PhyloZoo uses a hierarchical style system:

* **BaseStyle**: Common styling options for all plots
* **DMGraphStyle**: Extends BaseStyle for DirectedMultiGraph
* **MGraphStyle**: Extends BaseStyle for MixedMultiGraph
* **DNetStyle**: Extends DMGraphStyle for DirectedPhyNetwork
* **SDNetStyle**: Extends MGraphStyle for SemiDirectedPhyNetwork

BaseStyle Attributes
-------------------

All style classes inherit these common attributes from ``BaseStyle``:

* **node_color** (str, default='lightblue'): Color for nodes
* **node_size** (float, default=500.0): Size of nodes
* **edge_color** (str, default='gray'): Color for edges
* **edge_width** (float, default=2.0): Width of edges
* **with_labels** (bool, default=True): Whether to show node labels
* **label_offset** (float, default=0.1): Offset for labels from nodes
* **label_font_size** (float, default=10.0): Font size for labels
* **label_color** (str, default='black'): Color for labels

DNetStyle Attributes
--------------------

``DNetStyle`` extends ``DMGraphStyle`` with phylogenetic network-specific options:

* **leaf_color** (str, default='lightgreen'): Color for leaf nodes
* **hybrid_color** (str, default='salmon'): Color for hybrid nodes
* **leaf_size** (float, default=600.0): Size of leaf nodes
* **hybrid_edge_color** (str, default='red'): Color for hybrid edges

SDNetStyle Attributes
--------------------

``SDNetStyle`` extends ``MGraphStyle`` with semi-directed network-specific options:

* **leaf_color** (str, default='lightgreen'): Color for leaf nodes
* **hybrid_color** (str, default='salmon'): Color for hybrid nodes
* **leaf_size** (float, default=600.0): Size of leaf nodes
* **hybrid_edge_color** (str, default='red'): Color for hybrid edges

Using Styles
------------

**Get Default Style:**

.. code-block:: python

   from phylozoo.viz.dnetwork import default_style
   
   style = default_style()
   style.node_color = "blue"
   style.edge_color = "gray"

**Create Custom Style:**

.. code-block:: python

   from phylozoo.viz.dnetwork import DNetStyle
   from phylozoo.viz.sdnetwork import SDNetStyle
   
   # For directed networks
   d_style = DNetStyle(
       node_color="steelblue",
       leaf_color="lightgreen",
       hybrid_color="salmon",
       hybrid_edge_color="red",
       node_size=150,
       edge_width=2.0,
       label_font_size=10
   )
   
   # For semi-directed networks
   sd_style = SDNetStyle(
       node_color="lightblue",
       leaf_color="green",
       hybrid_color="orange",
       hybrid_edge_color="orange",
       node_size=150,
       edge_width=2.0
   )

**Modify Existing Style:**

.. code-block:: python

   from phylozoo.viz.dnetwork import default_style
   
   style = default_style()
   style.node_color = "blue"
   style.edge_color = "darkgray"
   style.node_size = 200
   style.edge_width = 3.0

**Use Style in Plotting:**

.. code-block:: python

   from phylozoo.viz import plot_dnetwork, DNetStyle
   
   style = DNetStyle(
       node_color="steelblue",
       hybrid_edge_color="red",
       node_size=150
   )
   
   plot_dnetwork(network, style=style, show=True)

Style Inheritance
-----------------

Styles follow an inheritance hierarchy:

* ``BaseStyle``: Common options (node_color, edge_color, etc.)
* ``DMGraphStyle(BaseStyle)``: For directed multigraphs
* ``MGraphStyle(BaseStyle)``: For mixed multigraphs
* ``DNetStyle(DMGraphStyle)``: For directed networks (adds leaf/hybrid options)
* ``SDNetStyle(MGraphStyle)``: For semi-directed networks (adds leaf/hybrid options)

This allows styles to be shared and extended appropriately for each graph/network type.

Color Options
------------

Colors can be specified as:

* **Named colors**: ``'red'``, ``'blue'``, ``'lightgreen'``, etc.
* **Hex codes**: ``'#FF0000'``, ``'#00FF00'``, etc.
* **RGB tuples**: ``(1.0, 0.0, 0.0)`` for red (normalized to 0-1)

Matplotlib color specifications are supported.

Node Types
----------

Different node types can have different colors and sizes:

* **Root nodes**: Use ``node_color`` and ``node_size``
* **Leaf nodes**: Use ``leaf_color`` and ``leaf_size``
* **Hybrid nodes**: Use ``hybrid_color`` and ``node_size``
* **Tree nodes**: Use ``node_color`` and ``node_size``

Edge Types
----------

Different edge types can have different colors:

* **Tree edges**: Use ``edge_color``
* **Hybrid edges**: Use ``hybrid_edge_color``

For DirectedPhyNetwork, all edges are directed and use ``edge_color`` or 
``hybrid_edge_color``. For SemiDirectedPhyNetwork, tree edges are undirected 
and use ``edge_color``, while hybrid edges are directed and use ``hybrid_edge_color``.

Complete Example
-----------------

.. code-block:: python

   from phylozoo import DirectedPhyNetwork
   from phylozoo.viz import plot_dnetwork, DNetStyle
   import matplotlib.pyplot as plt
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Create custom style
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
       label_color="black"
   )
   
   # Plot with custom style
   ax = plot_dnetwork(network, style=style, show=False)
   ax.figure.savefig("network.png", dpi=300, bbox_inches="tight")
   plt.close(ax.figure)

.. seealso::
   For layout options, see :doc:`Layouts <layouts>`. 
   For plotting functions, see :doc:`Plotting <plotting>`.
