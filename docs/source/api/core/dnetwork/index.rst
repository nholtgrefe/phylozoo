Directed Networks
=================

The :mod:`phylozoo.core.network.dnetwork` module provides the :class:`DirectedPhyNetwork` class
and related functions for working with directed phylogenetic networks. Directed phylogenetic
networks are rooted networks where all edges are directed. Internal nodes are either tree nodes
(in-degree 1, out-degree :math:`\geq 2`) or hybrid nodes (in-degree :math:`\geq 2`, out-degree 1).

Base Class
----------

.. automodule:: phylozoo.core.network.dnetwork.base
   :members:
   :show-inheritance:

Features
--------

.. automodule:: phylozoo.core.network.dnetwork.features
   :members:
   :show-inheritance:

Classifications
---------------

.. automodule:: phylozoo.core.network.dnetwork.classifications
   :members:
   :show-inheritance:

Transformations
---------------

.. automodule:: phylozoo.core.network.dnetwork.transformations
   :members:
   :show-inheritance:

Derivations
-----------

.. automodule:: phylozoo.core.network.dnetwork.derivations
   :members:
   :show-inheritance:

Conversions
-----------

.. automodule:: phylozoo.core.network.dnetwork.conversions
   :members:
   :show-inheritance:

Isomorphism
-----------

.. automodule:: phylozoo.core.network.dnetwork.isomorphism
   :members:
   :show-inheritance:

I/O
---

.. automodule:: phylozoo.core.network.dnetwork.io
   :members:
   :show-inheritance:

.. toctree::
   :maxdepth: 1
   
   generator
