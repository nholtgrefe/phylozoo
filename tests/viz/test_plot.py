"""
Tests for public plotting API in viz.

This test suite covers the plot_network function.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooLayoutError
from phylozoo.viz import plot_dnetwork
from phylozoo.viz.dnetwork import DNetStyle
from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID


class TestPlotNetwork:
    """Test plot_network function."""

    def test_simple_tree(self) -> None:
        """Test plotting a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        ax = plot_dnetwork(net, layout='pz-dag', trials=100)

        assert ax is not None

    def test_with_custom_style(self) -> None:
        """Test plotting with custom style."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        style = DNetStyle(node_color='blue', leaf_color='green')
        ax = plot_dnetwork(net, layout='pz-dag', style=style, trials=100)

        assert ax is not None

    def test_with_layout_kwargs(self) -> None:
        """Test plotting with layout parameters."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        ax = plot_dnetwork(net, layout='pz-dag', layer_gap=2.0, leaf_gap=3.0, trials=100)

        assert ax is not None

    def test_single_hybrid(self) -> None:
        """Test plotting network with hybrid node."""
        net = LEVEL_1_DNETWORK_SINGLE_HYBRID

        ax = plot_dnetwork(net, layout='pz-dag', trials=100)

        assert ax is not None

    def test_invalid_layout_raises(self) -> None:
        """Test that invalid layout raises PhyloZooLayoutError."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        with pytest.raises(PhyloZooLayoutError, match="not supported|Unknown|Unsupported"):
            plot_dnetwork(net, layout='invalid')
