dnetwork
--------

.. automodule:: phylozoo.viz.dnetwork
   :no-members:

Plotting
^^^^^^^^

.. autofunction:: phylozoo.viz.dnetwork.plot_dnetwork

Layouts
^^^^^^^

Each function returns a :class:`~phylozoo.viz.dnetwork.layout.base.DNetLayout`; the keyword
arguments are the layout parameters accepted by :func:`~phylozoo.viz.plot`.

.. autofunction:: phylozoo.viz.dnetwork.layout.cladogram.compute_pz_cladogram_layout

.. autofunction:: phylozoo.viz.dnetwork.layout.layered.compute_pz_layered_layout

.. autofunction:: phylozoo.viz.dnetwork.layout.radial.compute_pz_radial_layout

.. autofunction:: phylozoo.viz.dnetwork.layout.unrooted.compute_pz_unrooted_layout

.. autofunction:: phylozoo.viz.dnetwork.layout.nx.compute_nx_layout

.. autoclass:: phylozoo.viz.dnetwork.layout.base.DNetLayout
   :members: get_position, get_edge_route
   :exclude-members: network, positions, edge_routes, backbone_edges, reticulate_edges, algorithm, parameters
   :show-inheritance:

Styling
^^^^^^^

.. automodule:: phylozoo.viz.dnetwork.style
   :members:
   :show-inheritance:
