"""
Tests for layout computation in viz.

This test suite covers the DAG layout algorithm.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooLayoutError, PhyloZooValueError
from phylozoo.viz.dnetwork.layout import DAGLayout, compute_dag_layout
from tests.fixtures.directed_networks import (
    LEVEL_1_DNETWORK_SINGLE_HYBRID,
    LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
)


class TestDAGLayoutComputation:
    """Test DAG layout computation."""

    def test_simple_tree(self) -> None:
        """Test DAG layout on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        layout = compute_dag_layout(net, direction='TD', trials=100)

        assert isinstance(layout, DAGLayout)
        assert len(layout.positions) == 3
        assert len(layout.edge_routes) == 2
        assert layout.algorithm == 'dag'

    def test_single_hybrid(self) -> None:
        """Test DAG layout on network with single hybrid."""
        net = LEVEL_1_DNETWORK_SINGLE_HYBRID

        layout = compute_dag_layout(net, direction='TD', trials=100)

        assert isinstance(layout, DAGLayout)
        assert len(layout.positions) == net.number_of_nodes()
        assert len(layout.edge_routes) == net.number_of_edges()
        assert layout.algorithm == 'dag'

    def test_left_right_direction(self) -> None:
        """Test DAG layout with LR direction."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        layout = compute_dag_layout(net, direction='LR', trials=100)

        assert layout.algorithm == 'dag'
        assert len(layout.positions) == 3

    def test_empty_network_raises(self) -> None:
        """Test that empty network raises PhyloZooLayoutError."""
        net = DirectedPhyNetwork()

        with pytest.raises(PhyloZooLayoutError, match="empty network"):
            compute_dag_layout(net)

    def test_invalid_direction_raises(self) -> None:
        """Test that invalid direction raises PhyloZooValueError."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )

        with pytest.raises(PhyloZooValueError, match="direction must be"):
            compute_dag_layout(net, direction='invalid')
