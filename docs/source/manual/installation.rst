Installation
============

Installing PhyloZoo
-------------------

PhyloZoo can be installed using pip:

.. code-block:: bash

   pip install phylozoo

For development installation:

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
* Numba >= 0.56.0 (optional, for performance)

Optional Dependencies
---------------------

* **Visualization**: Matplotlib (for plotting networks)
* **Testing**: pytest, pytest-cov
* **Development**: mypy, black, ruff

Verifying Installation
----------------------

To verify that PhyloZoo is installed correctly:

.. code-block:: python

   >>> import phylozoo
   >>> print(phylozoo.__version__)
   0.1.0
