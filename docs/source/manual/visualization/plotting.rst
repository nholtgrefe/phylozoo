Plotting
========

The plotting system combines layout computation, styling, and rendering to create 
matplotlib figures. This page explains how the plotting process works.

Plotting Functions
------------------

PhyloZoo provides four main plotting functions:

* **plot_dnetwork()**: Plot DirectedPhyNetwork
* **plot_sdnetwork()**: Plot SemiDirectedPhyNetwork
* **plot_dmgraph()**: Plot DirectedMultiGraph
* **plot_mmgraph()**: Plot MixedMultiGraph

All plotting functions follow the same pattern:

1. **Layout Computation**: Compute node positions and edge routes
2. **Styling**: Apply style configuration
3. **Rendering**: Draw nodes, edges, and labels
4. **Display/Save**: Show or return the figure

Basic Usage
-----------

.. code-block:: python

   from phylozoo.viz import plot_dnetwork, plot_sdnetwork
   from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork
   
   # Plot with defaults
   dnet = DirectedPhyNetwork.load("network.enewick")
   plot_dnetwork(dnet, show=True)
   
   # Plot with custom layout
   sdnet = SemiDirectedPhyNetwork.load("network.enewick")
   plot_sdnetwork(sdnet, layout='spring', show=True)

Function Parameters
-------------------

All plotting functions accept these common parameters:

* **network/graph**: The network or graph to plot
* **layout** (str): Layout algorithm name (default varies by network type)
* **style** (Style | None): Styling configuration (uses default if None)
* **ax** (matplotlib.axes.Axes | None): Existing axes to plot on (creates new if None)
* **show** (bool): If True, automatically display the plot
* **\*\*layout_kwargs**: Additional parameters for layout computation

Layout Computation
-------------------

The plotting function first computes the layout:

1. **Check Layout Type**: Determine if it's a custom PhyloZoo layout (``pz-*``) or 
   a NetworkX/Graphviz layout
2. **Compute Positions**: Run the layout algorithm to get node positions
3. **Compute Routes**: Calculate edge routes (straight lines or curves)
4. **Create Layout Object**: Package positions and routes into a Layout object

The Layout object contains:

* **positions**: Dictionary mapping node IDs to (x, y) coordinates
* **edge_routes**: Dictionary mapping (u, v, key) to EdgeRoute objects
* **algorithm**: Name of the layout algorithm used
* **parameters**: Parameters used for layout computation

Rendering Process
-----------------

After layout computation, the plotting function renders the network:

1. **Render Edges**: Draw all edges with appropriate colors and styles
   * For DirectedPhyNetwork: All edges are drawn with arrows
   * For SemiDirectedPhyNetwork: Only hybrid edges (directed) have arrows
2. **Render Nodes**: Draw all nodes with appropriate colors and sizes
   * Different node types (root, leaf, hybrid, tree) use different colors
3. **Render Labels**: Add node labels if enabled
4. **Configure Axes**: Set aspect ratio and turn off axes

Edge Rendering
--------------

Edges are rendered based on their type:

* **Tree edges**: Use ``edge_color`` and ``edge_width``
* **Hybrid edges**: Use ``hybrid_edge_color`` and ``edge_width``
* **Parallel edges**: Offset slightly to avoid overlap

For DirectedPhyNetwork, all edges are directed and displayed with arrows pointing 
from parent to child.

For SemiDirectedPhyNetwork, only hybrid edges (which are directed) are displayed 
with arrows. Tree edges (which are undirected) are displayed as simple lines.

Node Rendering
--------------

Nodes are rendered based on their type:

* **Root nodes**: Use ``node_color`` and ``node_size``
* **Leaf nodes**: Use ``leaf_color`` and ``leaf_size``
* **Hybrid nodes**: Use ``hybrid_color`` and ``node_size``
* **Tree nodes**: Use ``node_color`` and ``node_size``

Node types are determined automatically from the network structure.

Label Rendering
---------------

If ``with_labels=True``, node labels are rendered:

* Labels are positioned with an offset from the node center
* Font size and color are controlled by style attributes
* For radial layouts, leaf labels may be positioned radially outward

Saving Figures
--------------

To save a plot to a file:

.. code-block:: python

   from phylozoo.viz import plot_dnetwork
   import matplotlib.pyplot as plt
   
   # Plot without showing
   ax = plot_dnetwork(network, show=False)
   
   # Save figure
   ax.figure.savefig("network.png", dpi=300, bbox_inches="tight")
   plt.close(ax.figure)

The function returns a matplotlib axes object, which provides access to the figure 
for saving or further customization.

Custom Axes
-----------

You can plot on an existing axes object:

.. code-block:: python

   import matplotlib.pyplot as plt
   from phylozoo.viz import plot_dnetwork
   
   # Create figure with subplots
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
   
   # Plot on first axes
   plot_dnetwork(network1, ax=ax1, show=False)
   
   # Plot on second axes
   plot_dnetwork(network2, ax=ax2, show=False)
   
   plt.tight_layout()
   plt.savefig("comparison.png", dpi=300)
   plt.close()

This allows you to create multi-panel figures or combine network plots with other 
visualizations.

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
       hybrid_edge_color="red",
       node_size=150
   )
   
   # Plot with custom layout and style
   ax = plot_dnetwork(
       network,
       layout='pz-dag',
       style=style,
       layer_gap=2.0,
       leaf_gap=1.5,
       show=False
   )
   
   # Save figure
   ax.figure.savefig("network.png", dpi=300, bbox_inches="tight")
   plt.close(ax.figure)

The plotting process is:

1. Layout ``'pz-dag'`` is computed with ``layer_gap=2.0`` and ``leaf_gap=1.5``
2. Custom style is applied
3. Network is rendered with the computed layout and style
4. Figure is saved to file

.. seealso::
   For layout algorithms, see :doc:`Layouts <layouts>`. 
   For styling options, see :doc:`Styling <styling>`. 
   For complete examples, see :doc:`Visualization Guide <viz>`.
