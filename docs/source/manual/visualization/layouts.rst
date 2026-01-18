Layout Algorithms
=================

Layout algorithms determine the positions of nodes and the routing of edges in the 
visualization. PhyloZoo provides both custom phylogenetic layouts and access to 
standard NetworkX and Graphviz layouts.

Custom PhyloZoo Layouts
------------------------

PhyloZoo includes specialized layout algorithms optimized for phylogenetic networks:

**For DirectedPhyNetwork:**

* **pz-dag** (default)
   Tree-backbone heuristic layout with crossing minimization. Optimizes child ordering 
   to minimize edge crossings. Suitable for all directed phylogenetic networks.
   
   Parameters:
   
   * ``layer_gap`` (float, default=1.5): Vertical spacing between layers
   * ``leaf_gap`` (float, default=1.0): Horizontal spacing between leaves
   * ``trials`` (int, default=2000): Number of trials for optimization
   * ``seed`` (int | None, default=None): Random seed for reproducibility
   * ``direction`` (str, default='TD'): Layout direction ('TD' for top-down, 'LR' for left-right)
   * ``x_scale`` (float, default=1.5): Horizontal scaling factor
   * ``y_scale`` (float, default=1.0): Vertical scaling factor
   
   .. code-block:: python
   
      from phylozoo.viz import plot_dnetwork
      
      plot_dnetwork(
          network,
          layout='pz-dag',
          layer_gap=2.0,
          leaf_gap=1.5,
          show=True
      )

**For SemiDirectedPhyNetwork:**

* **pz-radial**
   Radial (circular) layout for trees. Positions nodes in a circular arrangement with 
   the root at the center and leaves on the outer circle. Only works for trees 
   (``is_tree()`` must be True).
   
   Parameters:
   
   * ``radius`` (float, default=1.0): Outer radius for leaves
   * ``start_angle`` (float, default=0.0): Starting angle in radians
   * ``angle_direction`` (str, default='clockwise'): Direction of angle progression
   
   .. code-block:: python
   
      from phylozoo.viz import plot_sdnetwork
      
      plot_sdnetwork(
          network,
          layout='pz-radial',
          radius=2.0,
          show=True
      )

NetworkX Layouts
----------------

NetworkX provides a variety of general-purpose graph layout algorithms. These work 
for both DirectedPhyNetwork and SemiDirectedPhyNetwork:

* **spring**: Force-directed layout using Fruchterman-Reingold algorithm
* **circular**: Circular arrangement of nodes
* **kamada_kawai**: Force-directed layout using Kamada-Kawai algorithm
* **planar**: Planar layout (only for planar graphs)
* **random**: Random node positions
* **shell**: Concentric shell layout
* **spectral**: Spectral layout using eigenvectors
* **spiral**: Spiral arrangement
* **bipartite**: Bipartite layout (for bipartite graphs)

Parameters vary by algorithm. Common parameters include:

* ``k`` (float): Optimal distance between nodes (for spring layout)
* ``iterations`` (int): Number of iterations (for spring layout)
* ``seed`` (int | None): Random seed for reproducibility

.. code-block:: python

   from phylozoo.viz import plot_dnetwork, plot_sdnetwork
   
   # Use spring layout
   plot_dnetwork(network, layout='spring', k=2.0, iterations=50, show=True)
   
   # Use circular layout
   plot_sdnetwork(network, layout='circular', show=True)

Graphviz Layouts
----------------

Graphviz provides powerful graph layout algorithms. These require the ``pygraphviz`` 
package to be installed:

* **dot**: Hierarchical layout (default for directed graphs)
* **neato**: Spring-model layout
* **fdp**: Force-directed layout
* **sfdp**: Scalable force-directed layout
* **twopi**: Radial layout (default for SemiDirectedPhyNetwork)
* **circo**: Circular layout

.. code-block:: python

   from phylozoo.viz import plot_sdnetwork
   
   # Use twopi layout (default for SemiDirectedPhyNetwork)
   plot_sdnetwork(network, layout='twopi', show=True)
   
   # Use dot layout
   plot_sdnetwork(network, layout='dot', show=True)

.. note::
   Graphviz layouts require ``pygraphviz``. If not installed, NetworkX layouts will 
   be used instead. Install with: ``pip install pygraphviz``

Layout Selection
---------------

The choice of layout depends on your network structure and visualization goals:

* **For directed networks**: Use ``'pz-dag'`` for phylogenetic networks (default)
* **For semi-directed trees**: Use ``'pz-radial'`` for radial tree visualization
* **For general graphs**: Use NetworkX layouts like ``'spring'`` or ``'circular'``
* **For large networks**: Use Graphviz layouts like ``'twopi'`` or ``'sfdp'``

Layout Parameters
----------------

Layout-specific parameters can be passed directly to the plotting functions:

.. code-block:: python

   from phylozoo.viz import plot_dnetwork
   
   # Pass layout parameters
   plot_dnetwork(
       network,
       layout='pz-dag',
       layer_gap=2.0,
       leaf_gap=1.5,
       show=True
   )

All layout parameters are passed through ``**layout_kwargs`` to the layout computation 
function.

.. seealso::
   For more information on plotting functions, see :doc:`Plotting <plotting>`. 
   For styling options, see :doc:`Styling <styling>`.
