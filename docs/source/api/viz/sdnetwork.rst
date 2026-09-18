sdnetwork
---------

.. automodule:: phylozoo.viz.sdnetwork
   :no-members:

Plotting
^^^^^^^^

.. autofunction:: phylozoo.viz.sdnetwork.plot_sdnetwork

Layouts
^^^^^^^

Each function returns a :class:`~phylozoo.viz.sdnetwork.layout.base.SDNetLayout`; the keyword
arguments are the layout parameters accepted by :func:`~phylozoo.viz.plot`.

.. autofunction:: phylozoo.viz.sdnetwork.layout.unrooted.compute_pz_unrooted_layout

.. autofunction:: phylozoo.viz.sdnetwork.layout.radial.compute_pz_radial_layout

.. autofunction:: phylozoo.viz.sdnetwork.layout.nx.compute_nx_layout

.. autoclass:: phylozoo.viz.sdnetwork.layout.base.SDNetLayout
   :members: get_position, get_edge_route
   :exclude-members: network, positions, edge_routes, algorithm, parameters
   :show-inheritance:

Styling
^^^^^^^^

.. automodule:: phylozoo.viz.sdnetwork.style
   :members:
   :show-inheritance:
