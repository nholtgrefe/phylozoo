Overview
========

The :mod:`phylozoo.viz` module provides network and graph visualization. It consists of four submodules, each with styles and type-specific plotters:

- :mod:`phylozoo.viz.dnetwork` — :class:`~phylozoo.core.network.dnetwork.base.DirectedPhyNetwork`
- :mod:`phylozoo.viz.sdnetwork` — :class:`~phylozoo.core.network.sdnetwork.sd_phynetwork.SemiDirectedPhyNetwork`
- :mod:`phylozoo.viz.d_multigraph` — :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`
- :mod:`phylozoo.viz.m_multigraph` — :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`

There is a convenience function :func:`~phylozoo.viz.plot` that dispatches by object type, to one of the four submodules. It can be imported
directly from the :mod:`phylozoo.viz` module.

.. code-block:: python

   from phylozoo.viz import plot

A more in-depth explanation of plotting is available in the :doc:`Plotting <plotting>` documentation.
Styling options are explained in the :doc:`Styling <styling>` documentation.

See Also
--------

- :doc:`API reference <../../api/viz/index>` — Complete function signatures and detailed examples
