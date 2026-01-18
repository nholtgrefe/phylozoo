Overview
========

The visualization module (`phylozoo.viz`) provides a flexible plotting system for 
phylogenetic networks and graphs with support for multiple layout algorithms and 
customizable styling :cite:`PhyloZoo2024`.

The visualization system is organized into four main components:

* **Layout Algorithms**: Position nodes in 2D space
* **Styling**: Customize colors, sizes, and appearance
* **Plotting**: Render networks and graphs to matplotlib figures
* **Graph Support**: Plot underlying graph structures

Architecture
-----------

The visualization module uses a modular architecture:

1. **Layout Computation**: Layout algorithms compute node positions and edge routes
2. **Styling Configuration**: Style classes define visual appearance
3. **Rendering**: Plotting functions combine layouts and styles to create matplotlib figures

This separation allows for easy customization and extension of the visualization system.

Network Types
-------------

**DirectedPhyNetwork**
   Fully directed networks with all edges directed. Default layout: ``'pz-dag'``.
   All edges are displayed with arrows.

**SemiDirectedPhyNetwork**
   Networks with mixed directed and undirected edges. Default layout: ``'twopi'``.
   Only hybrid edges (directed edges) are displayed with arrows; undirected tree 
   edges are displayed without arrows.

Graph Types
-----------

**DirectedMultiGraph**
   Directed multigraphs supporting parallel edges. Uses standard NetworkX/Graphviz layouts.

**MixedMultiGraph**
   Mixed multigraphs with both directed and undirected edges. Uses standard NetworkX/Graphviz layouts.

Quick Start
-----------

.. code-block:: python

   from phylozoo.viz import plot_dnetwork, plot_sdnetwork
   from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork
   
   # Plot directed network (default layout: 'pz-dag')
   dnet = DirectedPhyNetwork.load("network.enewick")
   plot_dnetwork(dnet, show=True)
   
   # Plot semi-directed network (default layout: 'twopi')
   sdnet = SemiDirectedPhyNetwork.load("network.enewick")
   plot_sdnetwork(sdnet, show=True)

For more detailed information, see:

* :doc:`Layouts <layouts>`: Available layout algorithms and their parameters
* :doc:`Styling <styling>`: Customizing visual appearance
* :doc:`Plotting <plotting>`: How the plotting system works

.. seealso::
   For complete examples, see :doc:`Visualization Guide <viz>`. 
   For network analysis workflows, see :doc:`Network Analysis Workflow <../../tutorials/workflow_network_analysis>`.
