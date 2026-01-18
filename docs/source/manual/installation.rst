Installation
============

Installing PhyloZoo
-------------------

PhyloZoo can be installed using pip from PyPI:

.. code-block:: bash

   pip install phylozoo

For the latest development version, clone the repository and install in editable mode:

.. code-block:: bash

   git clone https://github.com/yourusername/phylozoo.git
   cd phylozoo
   pip install -e ".[dev]"

Requirements
------------

PhyloZoo requires:

* Python >= 3.10
* NumPy >= 1.20.0
* NetworkX >= 3.0.0
* Numba >= 0.56.0

Optional Dependencies
---------------------

* **Visualization**: Matplotlib >= 3.5.0 (for plotting networks and visualizations)
* **Testing**: pytest >= 7.0.0, pytest-cov >= 4.0.0
* **Development**: mypy >= 1.0.0, black >= 23.0.0, ruff >= 0.1.0
* **Documentation**: Sphinx >= 7.0.0, sphinx-autobuild, sphinx-rtd-theme, sphinxcontrib-napoleon

Platform Notes
--------------

* **Linux/macOS**: Full support for all features.
* **Windows**: Supported, but Numba JIT compilation may have issues on some configurations. If you encounter compilation errors, try installing from conda-forge or use the no-numba fallback (not recommended for performance).

Verifying Installation
----------------------

To verify that PhyloZoo is installed correctly:

.. code-block:: python

   >>> import phylozoo
   >>> print(phylozoo.__version__)
   0.1.0

   >>> from phylozoo import DirectedPhyNetwork
   >>> net = DirectedPhyNetwork(edges=[(1, 2), (2, 3)], nodes=[(3, {'label': 'A'})])
   >>> print(net.leaves)
   [3]

Troubleshooting
---------------

**Numba compilation errors**: If you see errors related to Numba JIT compilation, try:

.. code-block:: bash

   pip install --upgrade numba llvmlite

Or disable JIT compilation by setting the environment variable:

.. code-block:: bash

   export NUMBA_DISABLE_JIT=1

**Import errors**: Ensure you're using Python >= 3.10. Check that all dependencies are installed:

.. code-block:: bash

   pip check phylozoo

**Visualization not working**: Install Matplotlib:

.. code-block:: bash

   pip install matplotlib

Building Documentation
----------------------

To build the documentation locally:

.. code-block:: bash

   pip install -e ".[docs]"
   cd docs
   make html

Or with sphinx-autobuild for live reloading:

.. code-block:: bash

   sphinx-autobuild source build/html
