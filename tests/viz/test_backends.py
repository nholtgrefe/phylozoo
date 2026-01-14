"""
Tests for backend system in viz.

This test suite covers the backend registry and matplotlib backend.
"""

import pytest

from phylozoo.utils.exceptions import PhyloZooBackendError
from phylozoo.viz.dnetwork.backends import (
    Backend,
    get_backend,
    register_backend,
)


class TestBackendRegistry:
    """Test backend registry system."""

    def test_register_and_get_backend(self) -> None:
        """Test registering and retrieving a backend."""
        # Create a dummy backend class
        class DummyBackend(Backend):
            def create_figure(self, **kwargs):
                return None

            def create_axes(self, figure=None, **kwargs):
                return None

            def render_edge(self, route, style, layout):
                return {}

            def render_node(self, node, position, node_type, style):
                return None

            def add_label(self, position, text, style):
                return None

            def show(self):
                pass

            def save(self, path, **kwargs):
                pass

        # Register the backend
        register_backend('dummy', DummyBackend)

        # Retrieve it
        BackendClass = get_backend('dummy')
        assert BackendClass == DummyBackend

    def test_get_nonexistent_backend_raises(self) -> None:
        """Test that getting nonexistent backend raises PhyloZooBackendError."""
        with pytest.raises(PhyloZooBackendError, match="not registered"):
            get_backend('nonexistent_backend')

    def test_matplotlib_backend_registered(self) -> None:
        """Test that matplotlib backend is registered by default."""
        # Import should register it
        from phylozoo.viz.dnetwork.backends import matplotlib  # noqa: F401

        BackendClass = get_backend('matplotlib')
        assert BackendClass is not None


class TestMatplotlibBackend:
    """Test matplotlib backend implementation."""

    def test_create_figure(self) -> None:
        """Test figure creation."""
        from phylozoo.viz.dnetwork.backends import get_backend

        BackendClass = get_backend('matplotlib')
        backend = BackendClass()

        fig = backend.create_figure()
        assert fig is not None
        assert backend.figure is not None

    def test_create_axes(self) -> None:
        """Test axes creation."""
        from phylozoo.viz.dnetwork.backends import get_backend

        BackendClass = get_backend('matplotlib')
        backend = BackendClass()

        backend.create_figure()
        ax = backend.create_axes()
        assert ax is not None
        assert backend.axes is not None

