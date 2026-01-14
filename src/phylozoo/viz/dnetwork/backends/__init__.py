"""Backend implementations for DirectedPhyNetwork plotting."""

from .base import Backend, get_backend, register_backend

# Import backends to register them
from . import matplotlib  # noqa: F401

# Try to import pyqtgraph backend (optional)
try:
    from . import pyqtgraph  # noqa: F401
except ImportError:
    pass

__all__ = ['Backend', 'get_backend', 'register_backend']

