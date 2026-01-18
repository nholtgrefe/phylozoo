Visualization Module
====================

The visualization module (`phylozoo.viz`) provides a flexible plotting system for 
phylogenetic networks and graphs with support for multiple layout algorithms and 
customizable styling :cite:`PhyloZoo2024`.

Basic Plotting
--------------

The simplest way to plot a network is using the plotting functions:

.. code-block:: python

   from phylozoo.viz import plot_dnetwork, plot_sdnetwork
   from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork
   
   # Plot directed network (default layout: 'pz-dag')
   dnet = DirectedPhyNetwork.load("network.enewick")
   plot_dnetwork(dnet, show=True)
   
   # Plot semi-directed network (default layout: 'twopi')
   sdnet = SemiDirectedPhyNetwork.load("network.enewick")
   plot_sdnetwork(sdnet, show=True)

The ``show=True`` parameter displays the plot immediately. Set ``show=False`` to 
return a figure object for further customization or saving.

Layout Algorithms
-----------------

PhyloZoo provides multiple layout algorithms for positioning nodes:

**For DirectedPhyNetwork (default: 'pz-dag'):**

.. code-block:: python

   from phylozoo.viz import plot_dnetwork
   
   # Use default pz-dag layout
   plot_dnetwork(network, show=True)
   
   # Use NetworkX or Graphviz layouts
   plot_dnetwork(network, layout='spring', show=True)
   plot_dnetwork(network, layout='dot', show=True)

**For SemiDirectedPhyNetwork (default: 'twopi'):**

.. code-block:: python

   from phylozoo.viz import plot_sdnetwork
   
   # Use default twopi layout
   plot_sdnetwork(network, show=True)
   
   # Use custom PhyloZoo layouts
   plot_sdnetwork(network, layout='pz-radial', show=True)
   plot_sdnetwork(network, layout='pz-planar', show=True)
   
   # Use NetworkX layouts
   plot_sdnetwork(network, layout='spring', show=True)

Layout parameters can be passed directly to the plotting functions:

.. code-block:: python

   # Pass layout parameters directly
   plot_dnetwork(
       network,
       node_gap=1.0,
       layer_gap=1.5,
       show=True
   )

Styling
-------

Networks can be styled using style classes:

.. code-block:: python

   from phylozoo.viz import plot_dnetwork, plot_sdnetwork, DNetStyle, SDNetStyle, default_style
   
   # Get default style
   style = default_style()
   
   # Modify style
   style.node_color = "blue"
   style.edge_color = "gray"
   style.node_size = 100
   
   # Plot with custom style
   plot_dnetwork(network, style=style, show=True)

Or create a new style:

.. code-block:: python

   # Create custom style for directed networks
   style = DNetStyle(
       node_color="blue",
       edge_color="gray",
       node_size=100,
       edge_width=2.0
   )
   
   plot_dnetwork(network, style=style, show=True)
   
   # Create custom style for semi-directed networks
   sd_style = SDNetStyle(
       node_color="green",
       edge_color="darkgray",
       hybrid_edge_color="red",
       node_size=100
   )
   
   plot_sdnetwork(sdnet, style=sd_style, show=True)

Edge Rendering
--------------

PhyloZoo automatically handles edge directionality:

* **DirectedPhyNetwork**: All edges are directed and displayed with arrows
* **SemiDirectedPhyNetwork**: Only hybrid edges (directed edges) are displayed with arrows; 
  undirected tree edges are displayed without arrows

Graph Plotting
--------------

The visualization module also provides functions for plotting the underlying graph 
structures:

.. code-block:: python

   from phylozoo.viz.graphs import plot_dmgraph, plot_mmgraph
   from phylozoo.core.primitives.dmultigraph import DirectedMultiGraph
   
   # Plot directed multigraph
   dmg = DirectedMultiGraph(edges=[(1, 2), (1, 3)])
   plot_dmgraph(dmg, show=True)

Saving Figures
-------------

Plots can be saved to files:

.. code-block:: python

   from phylozoo.viz.dnetwork import plot_network
   import matplotlib.pyplot as plt
   
   # Plot and save
   fig = plot_network(network, show=False)
   fig.savefig("network.png", dpi=300, bbox_inches="tight")
   plt.close(fig)

Complete Example
----------------

Here's a complete example with custom styling and layout:

.. code-block:: python

   from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork
   from phylozoo.viz import plot_dnetwork, plot_sdnetwork, DNetStyle, SDNetStyle
   import matplotlib.pyplot as plt
   
   # Plot directed network
   dnet = DirectedPhyNetwork.load("network.enewick")
   
   # Create custom style
   style = DNetStyle(
       node_color="steelblue",
       edge_color="darkgray",
       hybrid_edge_color="red",
       node_size=150,
       edge_width=2.0,
       label_fontsize=10
   )
   
   # Plot with custom layout and style
   ax = plot_dnetwork(
       dnet,
       layout='pz-dag',
       style=style,
       node_gap=1.2,
       layer_gap=1.5,
       show=False
   )
   
   # Save
   ax.figure.savefig("dnetwork.png", dpi=300, bbox_inches="tight")
   plt.close(ax.figure)
   
   # Plot semi-directed network
   sdnet = SemiDirectedPhyNetwork.load("network.enewick")
   
   # Create custom style
   sd_style = SDNetStyle(
       node_color="green",
       edge_color="darkgray",
       hybrid_edge_color="orange",
       node_size=150,
       edge_width=2.0
   )
   
   # Plot with twopi layout (default)
   ax2 = plot_sdnetwork(sdnet, style=sd_style, show=False)
   ax2.figure.savefig("sdnetwork.png", dpi=300, bbox_inches="tight")
   plt.close(ax2.figure)

Available Functions
-------------------

**Network Plotting:**

* **plot_dnetwork(network, layout='pz-dag', style=None, ax=None, show=False, **layout_kwargs)** - 
  Plot a DirectedPhyNetwork. Main plotting function with customizable layout and styling. 
  Default layout is 'pz-dag'. All edges are displayed with arrows. Returns matplotlib axes object. 
  See :func:`phylozoo.viz.dnetwork.plot_dnetwork`.

* **plot_sdnetwork(network, layout='twopi', style=None, ax=None, show=False, **layout_kwargs)** - 
  Plot a SemiDirectedPhyNetwork. Main plotting function with customizable layout and styling. 
  Default layout is 'twopi'. Only hybrid edges (directed edges) are displayed with arrows; 
  undirected tree edges are displayed without arrows. Returns matplotlib axes object. 
  See :func:`phylozoo.viz.sdnetwork.plot_sdnetwork`.

**Layout Algorithms:**

* **Custom PhyloZoo layouts for DirectedPhyNetwork:**
  - 'pz-dag': DAG-based layout (default)

* **Custom PhyloZoo layouts for SemiDirectedPhyNetwork:**
  - 'pz-radial': Radial layout for trees
  - 'pz-planar': Planar circular layout for planar networks

* **NetworkX layouts:** 'spring', 'circular', 'kamada_kawai', 'planar', 'random', etc.

* **Graphviz layouts:** 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo' (requires pygraphviz)

**Styling:**

* **DNetStyle** - Style configuration class for DirectedPhyNetwork. Attributes include 
  node_color, edge_color, hybrid_edge_color, node_size, edge_width, label_fontsize, etc. 
  See :class:`phylozoo.viz.dnetwork.DNetStyle`.

* **SDNetStyle** - Style configuration class for SemiDirectedPhyNetwork. Attributes include 
  node_color, edge_color, hybrid_edge_color, node_size, edge_width, label_fontsize, etc. 
  See :class:`phylozoo.viz.sdnetwork.SDNetStyle`.

* **default_style()** - Get default network style. Returns DNetStyle instance with 
  default values. See :func:`phylozoo.viz.dnetwork.default_style`.

**Graph Plotting:**

* **plot_dmgraph(graph, layout='spring', style=None, ax=None, show=False, **layout_kwargs)** - 
  Plot DirectedMultiGraph. Returns matplotlib axes object. 
  See :func:`phylozoo.viz.graphs.dmgraph.plot_dmgraph`.

* **plot_mmgraph(graph, layout='spring', style=None, ax=None, show=False, **layout_kwargs)** - 
  Plot MixedMultiGraph. Returns matplotlib axes object. 
  See :func:`phylozoo.viz.graphs.mmgraph.plot_mmgraph`.

.. note::
   The visualization module uses a modular architecture with separate layout, styling, 
   and rendering components. This allows for easy customization and extension.

.. tip::
   Use ``show=False`` to get the axes object for further customization before displaying 
   or saving. The default layouts ('pz-dag' for DirectedPhyNetwork and 'twopi' for 
   SemiDirectedPhyNetwork) are optimized for phylogenetic network visualization.

.. warning::
   Graphviz layouts require pygraphviz installation. If not available, NetworkX layouts 
   will be used instead.

.. seealso::
   For network analysis workflows that include visualization, see 
   :doc:`Network Analysis Workflow <../../tutorials/workflow_network_analysis>`. For I/O operations, 
   see :doc:`I/O <io>`.
