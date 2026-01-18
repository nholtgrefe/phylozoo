Installation
============

Installing PhyloZoo
-------------------

PhyloZoo can be installed using pip from PyPI:

.. code-block:: bash

   pip install phylozoo

For the latest development version, clone the repository and install in editable mode:

.. code-block:: bash

   git clone https://github.com/nholtgrefe/phylozoo.git
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

Verifying Installation
----------------------

To verify that PhyloZoo is installed correctly:

.. code-block:: python

   >>> import phylozoo
   >>> print(phylozoo.__version__)
   0.1.0

   >>> from phylozoo import DirectedPhyNetwork
   >>> net = DirectedPhyNetwork(edges=[(1, 2), (1, 3)], nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})])
   >>> print(net.leaves)
   [2, 3]

Troubleshooting
---------------

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
