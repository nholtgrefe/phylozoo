Visualization Module
====================

The visualization module (`phylozoo.viz`) provides a flexible plotting system for 
phylogenetic networks and graphs with support for multiple backends and customizable 
styling :cite:`PhyloZoo2024`.

Basic Plotting
--------------

The simplest way to plot a network is using the ``plot_network()`` function:

.. code-block:: python

   from phylozoo.viz.dnetwork import plot_network
   from phylozoo import DirectedPhyNetwork
   
   # Load and plot network
   network = DirectedPhyNetwork.load("network.enewick")
   plot_network(network, show=True)

The ``show=True`` parameter displays the plot immediately. Set ``show=False`` to 
return a figure object for further customization or saving.

Plotting Trees
--------------

For tree networks, you can use the specialized ``plot_tree()`` function:

.. code-block:: python

   from phylozoo.viz.dnetwork import plot_tree
   
   # Plot tree (optimized for tree structures)
   plot_tree(network, show=True)

Layout Algorithms
-----------------

PhyloZoo provides layout algorithms for positioning nodes:

.. code-block:: python

   from phylozoo.viz.dnetwork import compute_dag_layout, plot_network
   
   # Compute custom layout
   layout = compute_dag_layout(
       network,
       node_gap=1.0,
       layer_gap=1.5
   )
   
   # Plot with custom layout
   plot_network(network, layout=layout, show=True)

Layout parameters can be passed directly to ``plot_network()``:

.. code-block:: python

   # Pass layout parameters directly
   plot_network(
       network,
       node_gap=1.0,
       layer_gap=1.5,
       show=True
   )

Styling
-------

Networks can be styled using the ``NetworkStyle`` class:

.. code-block:: python

   from phylozoo.viz.dnetwork import NetworkStyle, plot_network, default_style
   
   # Get default style
   style = default_style()
   
   # Modify style
   style.node_color = "blue"
   style.edge_color = "gray"
   style.node_size = 100
   
   # Plot with custom style
   plot_network(network, style=style, show=True)

Or create a new style:

.. code-block:: python

   # Create custom style
   style = NetworkStyle(
       node_color="blue",
       edge_color="gray",
       node_size=100,
       edge_width=2.0
   )
   
   plot_network(network, style=style, show=True)

Backends
--------

PhyloZoo supports multiple plotting backends:

* **Matplotlib**: Default backend, produces static images (always available)
* **PyQtGraph**: Optional backend for interactive plots (requires installation)

The backend is automatically selected based on availability. Specify backend explicitly:

.. code-block:: python

   # Use matplotlib backend (default)
   plot_network(network, backend="matplotlib", show=True)
   
   # Use PyQtGraph backend (if available)
   plot_network(network, backend="pyqtgraph", show=True)

Graph Plotting
--------------

The visualization module also provides functions for plotting the underlying graph 
structures:

.. code-block:: python

   from phylozoo.viz.graphs import plot_directed_multigraph, plot_mixed_multigraph
   from phylozoo.core.primitives.dmultigraph import DirectedMultiGraph
   
   # Plot directed multigraph
   dmg = DirectedMultiGraph(edges=[(1, 2), (1, 3)])
   plot_directed_multigraph(dmg, show=True)

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

   from phylozoo import DirectedPhyNetwork
   from phylozoo.viz.dnetwork import plot_network, NetworkStyle
   
   # Load network
   network = DirectedPhyNetwork.load("network.enewick")
   
   # Create custom style
   style = NetworkStyle(
       node_color="steelblue",
       edge_color="darkgray",
       node_size=150,
       edge_width=2.0,
       label_fontsize=10
   )
   
   # Plot with custom layout and style
   fig = plot_network(
       network,
       style=style,
       node_gap=1.2,
       layer_gap=1.5,
       show=False
   )
   
   # Save
   fig.savefig("network.png", dpi=300, bbox_inches="tight")

Available Functions
-------------------

**Network Plotting:**

* **plot_network(network, layout='dag', style=None, backend='matplotlib', ax=None, show=False, **layout_kwargs)** - 
  Plot a DirectedPhyNetwork. Main plotting function with customizable layout, styling, 
  and backend. Returns backend-specific figure/axes object. See 
  :func:`phylozoo.viz.dnetwork.plot_network`.

* **plot_tree(network, layout='dag', style=None, backend='matplotlib', ax=None, show=False, **layout_kwargs)** - 
  Plot a tree network. Optimized for tree structures. Same parameters as plot_network. 
  See :func:`phylozoo.viz.dnetwork.plot_tree`.

**Layout:**

* **compute_dag_layout(network, node_gap=1.0, layer_gap=1.0, iterations=50, curve_strength=0.5)** - 
  Compute DAG layout for network. Positions nodes in layers with automatic spacing. 
  Returns layout object. See :func:`phylozoo.viz.dnetwork.compute_dag_layout`.

**Styling:**

* **NetworkStyle** - Style configuration class. Attributes include node_color, edge_color, 
  node_size, edge_width, label_fontsize, etc. See :class:`phylozoo.viz.dnetwork.NetworkStyle`.

* **default_style()** - Get default network style. Returns NetworkStyle instance with 
  default values. See :func:`phylozoo.viz.dnetwork.default_style`.

**Backends:**

* **get_backend(name)** - Get backend class by name. Returns backend class. See 
  :func:`phylozoo.viz.dnetwork.get_backend`.

* **register_backend(name, backend_class)** - Register custom backend. Allows extending 
  visualization with custom backends. See :func:`phylozoo.viz.dnetwork.register_backend`.

**Graph Plotting:**

* **plot_directed_multigraph(graph, show=False)** - Plot DirectedMultiGraph. Returns 
  matplotlib figure. See :func:`phylozoo.viz.graphs.plot_directed_multigraph`.

* **plot_mixed_multigraph(graph, show=False)** - Plot MixedMultiGraph. Returns matplotlib 
  figure. See :func:`phylozoo.viz.graphs.plot_mixed_multigraph`.

.. note::
   The visualization module uses a modular architecture with separate layout, styling, 
   and rendering components. This allows for easy customization and extension.

.. tip::
   Use ``show=False`` to get the figure object for further customization before displaying 
   or saving.

.. warning::
   PyQtGraph backend requires additional installation. If not available, matplotlib will 
   be used automatically.

.. seealso::
   For network analysis workflows that include visualization, see 
   :doc:`Network Analysis Workflow <workflow_network_analysis>`. For I/O operations, 
   see :doc:`I/O <io>`.
