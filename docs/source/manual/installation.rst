Installation
============

Installing PhyloZoo
-------------------

PhyloZoo is a Python package that runs on `Python <https://www.python.org/>`_ (>= 3.10).
Install PhyloZoo using ``pip`` from `PyPI <https://pypi.org/project/phylozoo/>`_. Choose one of:

* **Minimal** — Core only: NumPy, Numba, NetworkX.

  .. code-block:: bash

     pip install phylozoo

* **With plotting (recommended)** — Adds Matplotlib for network visualization.

  .. code-block:: bash

     pip install phylozoo[viz]

For development and contributing to PhyloZoo, install the latest source version in
editable mode:

.. code-block:: bash

   git clone https://github.com/nholtgrefe/phylozoo.git
   cd phylozoo
   pip install -e ".[dev]"

Requirements
^^^^^^^^^^^^

PhyloZoo keeps its core dependencies minimal. The mandatory requirements are:

* `NumPy <https://numpy.org/>`_ >= 1.20.0 (for numerical operations)
* `NetworkX <https://networkx.org/>`_ >= 3.0.0 (for graph operations)
* `Numba <https://numba.pydata.org/>`_ >= 0.56.0 (for JIT compilation of computationally intensive algorithms)

Optional (install via extras):

* `Matplotlib <https://matplotlib.org/>`_ >= 3.5.0 (for plotting; use ``phylozoo[viz]``)

Graphviz layouts (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

PhyloZoo's own layouts (``pz-cladogram``, ``pz-layered``, ``pz-radial``, ``pz-unrooted``) are the
recommended way to draw networks and need nothing beyond ``phylozoo[viz]``. The Graphviz layouts
(``dot``, ``neato``, ``twopi``, ...) can additionally be used if `PyGraphviz <https://pygraphviz.github.io/>`_
is installed; it needs the Graphviz system library (e.g. ``apt install graphviz graphviz-dev`` on
Debian/Ubuntu, ``brew install graphviz`` on macOS) and then

.. code-block:: bash

   pip install pygraphviz

See the `PyGraphviz installation guide <https://pygraphviz.github.io/documentation/stable/install.html>`_
for platform details.

Verifying Installation
----------------------

To verify that PhyloZoo is installed correctly, you can import it and print the version.
The latest version is |version|.

.. code-block:: python

   >>> import phylozoo
   >>> print(phylozoo.__version__)
   x.y.z  # your installed version

   >>> from phylozoo import DirectedPhyNetwork
   >>> net = DirectedPhyNetwork(edges=[(1, 2), (1, 3)], nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})])
   >>> print(net.leaves)
   [2, 3]


Building Documentation
----------------------

To build the documentation locally, install the optional documentation dependencies using
the ``docs`` extra:

.. code-block:: bash

   pip install -e ".[docs]"
   cd docs
   make html

This extra installs the main documentation tools:

* `Sphinx <https://www.sphinx-doc.org/>`_ >= 7.0.0
* `sphinx-autobuild <https://github.com/executablebooks/sphinx-autobuild>`_
* `sphinx-rtd-theme <https://sphinx-rtd-theme.readthedocs.io/>`_
* `sphinxcontrib-napoleon <https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html>`_

Or with ``sphinx-autobuild`` for live reloading:

.. code-block:: bash

   sphinx-autobuild source build/html


Troubleshooting
---------------

**Import errors**: Ensure you're using Python >= 3.10. Check that all dependencies are installed:

.. code-block:: bash

   pip check phylozoo

**Visualization not working / ``PhyloZooImportError``**: The viz module requires Matplotlib.
Install the viz extra:

.. code-block:: bash

   pip install phylozoo[viz]

**Graphviz layouts (dot, neato, fdp, etc.) not working**: These need the Graphviz system library
and the ``pygraphviz`` package, see *Graphviz layouts (optional)* above. The PhyloZoo layouts do
not need them.
