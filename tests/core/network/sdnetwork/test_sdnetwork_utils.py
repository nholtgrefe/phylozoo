"""
Tests for semi-directed network utility functions.

This test suite covers internal utility functions for network transformations,
including edge attribute splitting for edge subdivision.
"""

import pytest

from phylozoo.core.network.sdnetwork._utils import (
    _split_attrs_for_subdividing_edge,
    _subdivide_edge,
)


class TestSplitAttrsForSubdividingEdge:
    """Test cases for _split_attrs_for_subdividing_edge function."""

    def test_basic_with_branch_length_and_gamma(self) -> None:
        """Test splitting edge with branch_length and gamma."""
        edge_data = {'branch_length': 1.0, 'gamma': 0.6, 'bootstrap': 0.95}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.5)
        
        # First edge should have branch_length only (split at 0.5)
        assert first_attrs == {'branch_length': 0.5}
        
        # Second edge should have branch_length and gamma
        assert second_attrs == {'branch_length': 0.5, 'gamma': 0.6}
        
        # Bootstrap should be removed
        assert 'bootstrap' not in first_attrs
        assert 'bootstrap' not in second_attrs

    def test_different_subdivision_location(self) -> None:
        """Test splitting with different subdivision location."""
        edge_data = {'branch_length': 2.0}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.3)
        
        # First edge: 2.0 * 0.3 = 0.6
        assert first_attrs == {'branch_length': 0.6}
        
        # Second edge: 2.0 * 0.7 = 1.4
        assert second_attrs == {'branch_length': 1.4}

    def test_subdivision_at_edges(self) -> None:
        """Test splitting at edge cases (0.0 and 1.0)."""
        edge_data = {'branch_length': 1.0}
        
        # At 0.0: first edge gets nothing, second gets all
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.0)
        assert first_attrs == {'branch_length': 0.0}
        assert second_attrs == {'branch_length': 1.0}
        
        # At 1.0: first edge gets all, second gets nothing
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 1.0)
        assert first_attrs == {'branch_length': 1.0}
        assert second_attrs == {'branch_length': 0.0}

    def test_only_gamma_no_branch_length(self) -> None:
        """Test splitting edge with only gamma, no branch_length."""
        edge_data = {'gamma': 0.4, 'bootstrap': 0.9}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data)
        
        # First edge gets nothing
        assert first_attrs == {}
        
        # Second edge gets gamma only
        assert second_attrs == {'gamma': 0.4}
        
        # Bootstrap should be removed
        assert 'bootstrap' not in first_attrs
        assert 'bootstrap' not in second_attrs

    def test_only_branch_length_no_gamma(self) -> None:
        """Test splitting edge with only branch_length, no gamma."""
        edge_data = {'branch_length': 1.5}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.5)
        
        assert first_attrs == {'branch_length': 0.75}
        assert second_attrs == {'branch_length': 0.75}
        assert 'gamma' not in first_attrs
        assert 'gamma' not in second_attrs

    def test_empty_edge_data(self) -> None:
        """Test splitting empty edge data."""
        edge_data = {}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data)
        
        assert first_attrs == {}
        assert second_attrs == {}

    def test_none_edge_data(self) -> None:
        """Test splitting None edge data (treated as empty)."""
        edge_data = None
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data)
        
        assert first_attrs == {}
        assert second_attrs == {}

    def test_other_attributes_removed(self) -> None:
        """Test that other attributes (bootstrap, etc.) are removed."""
        edge_data = {
            'branch_length': 1.0,
            'gamma': 0.5,
            'bootstrap': 0.95,
            'custom_attr': 'value',
            'weight': 10.0
        }
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data)
        
        # Only branch_length and gamma should remain
        assert set(first_attrs.keys()) == {'branch_length'}
        assert set(second_attrs.keys()) == {'branch_length', 'gamma'}
        
        # Other attributes should be removed
        assert 'bootstrap' not in first_attrs
        assert 'bootstrap' not in second_attrs
        assert 'custom_attr' not in first_attrs
        assert 'custom_attr' not in second_attrs
        assert 'weight' not in first_attrs
        assert 'weight' not in second_attrs

    def test_invalid_subdivision_location_below_zero(self) -> None:
        """Test that subdivision_location below 0.0 raises ValueError."""
        edge_data = {'branch_length': 1.0}
        
        with pytest.raises(ValueError, match="subdivision_location must be in \\[0.0, 1.0\\]"):
            _split_attrs_for_subdividing_edge(edge_data, -0.1)

    def test_invalid_subdivision_location_above_one(self) -> None:
        """Test that subdivision_location above 1.0 raises ValueError."""
        edge_data = {'branch_length': 1.0}
        
        with pytest.raises(ValueError, match="subdivision_location must be in \\[0.0, 1.0\\]"):
            _split_attrs_for_subdividing_edge(edge_data, 1.5)

    def test_default_subdivision_location(self) -> None:
        """Test that default subdivision_location is 0.5."""
        edge_data = {'branch_length': 1.0}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data)
        
        # Should split evenly at 0.5
        assert first_attrs == {'branch_length': 0.5}
        assert second_attrs == {'branch_length': 0.5}

    def test_gamma_only_on_second_edge(self) -> None:
        """Test that gamma always goes to second edge, never first."""
        edge_data = {'branch_length': 1.0, 'gamma': 0.7}
        
        # Test with different subdivision locations
        for loc in [0.0, 0.1, 0.5, 0.9, 1.0]:
            first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, loc)
            assert 'gamma' not in first_attrs
            assert 'gamma' in second_attrs
            assert second_attrs['gamma'] == 0.7

    def test_branch_length_precision(self) -> None:
        """Test that branch_length splitting maintains precision."""
        edge_data = {'branch_length': 1.0}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.333333)
        
        # Should handle floating point precision correctly
        assert abs(first_attrs['branch_length'] - 0.333333) < 1e-6
        assert abs(second_attrs['branch_length'] - 0.666667) < 1e-6
        # Sum should equal original
        assert abs(first_attrs['branch_length'] + second_attrs['branch_length'] - 1.0) < 1e-6

    def test_zero_branch_length(self) -> None:
        """Test splitting edge with zero branch_length."""
        edge_data = {'branch_length': 0.0, 'gamma': 0.5}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.5)
        
        assert first_attrs == {'branch_length': 0.0}
        assert second_attrs == {'branch_length': 0.0, 'gamma': 0.5}

    def test_negative_branch_length(self) -> None:
        """Test splitting edge with negative branch_length (edge case)."""
        edge_data = {'branch_length': -1.0}
        first_attrs, second_attrs = _split_attrs_for_subdividing_edge(edge_data, 0.5)
        
        # Function should still work (doesn't validate branch_length sign)
        assert first_attrs == {'branch_length': -0.5}
        assert second_attrs == {'branch_length': -0.5}


class TestSubdivideEdge:
    """Test cases for _subdivide_edge function."""

    def test_subdivide_undirected_edge(self) -> None:
        """Test subdividing an undirected edge."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 3, 'v': 1, 'key': 0, 'branch_length': 1.0},
                (3, 2), (3, 4)  # Node 3 needs degree >= 3
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        graph, subdiv_node = _subdivide_edge(net, 3, 1, 0, 0.5)
        
        # Check subdivision node exists
        assert subdiv_node in graph.nodes()
        
        # Check original edge is removed
        assert not graph.has_edge(3, 1, key=0)
        
        # Check new edges exist (undirected)
        assert graph.has_edge(3, subdiv_node)
        assert graph.has_edge(subdiv_node, 1)
        # Verify they are undirected
        assert graph._undirected.has_edge(3, subdiv_node)
        assert graph._undirected.has_edge(subdiv_node, 1)
        
        # Check attributes are split
        edge_data_1 = graph._undirected[3][subdiv_node][0]
        edge_data_2 = graph._undirected[subdiv_node][1][0]
        assert edge_data_1.get('branch_length') == 0.5
        assert edge_data_2.get('branch_length') == 0.5

    def test_subdivide_directed_edge(self) -> None:
        """Test subdividing a directed edge."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.utils.validation import no_validation
        
        # Use a valid network structure similar to test fixtures
        with no_validation():
            net = SemiDirectedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 6, 'v': 4, 'gamma': 0.4},
                    {'u': 3, 'v': 1, 'key': 0, 'branch_length': 1.0, 'gamma': 0.6}
                ],
                undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11), (3, 1)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        
        graph, subdiv_node = _subdivide_edge(net, 3, 1, 0, 0.5)
        
        # Check subdivision node exists
        assert subdiv_node in graph.nodes()
        
        # Check original edge is removed
        assert not graph._directed.has_edge(3, 1, key=0)
        
        # Check new edges exist: u-w (undirected) and w->v (directed)
        assert graph._undirected.has_edge(3, subdiv_node)
        assert graph._directed.has_edge(subdiv_node, 1)
        
        # Check attributes are split correctly
        edge_data_1 = graph._undirected[3][subdiv_node][0]
        edge_data_2 = graph._directed[subdiv_node][1][0]
        assert edge_data_1.get('branch_length') == 0.5
        assert edge_data_2.get('branch_length') == 0.5
        # Gamma should be on second edge only
        assert 'gamma' not in edge_data_1
        assert edge_data_2.get('gamma') == 0.6

    def test_subdivide_edge_different_location(self) -> None:
        """Test subdividing edge with different subdivision location."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 3, 'v': 1, 'key': 0, 'branch_length': 2.0},
                (3, 2), (3, 4)  # Node 3 needs degree >= 3
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        graph, subdiv_node = _subdivide_edge(net, 3, 1, 0, 0.3)
        
        # Check attributes are split at 0.3
        edge_data_1 = graph._undirected[3][subdiv_node][0]
        edge_data_2 = graph._undirected[subdiv_node][1][0]
        assert edge_data_1.get('branch_length') == 0.6  # 2.0 * 0.3
        assert edge_data_2.get('branch_length') == 1.4  # 2.0 * 0.7

    def test_subdivide_edge_not_found(self) -> None:
        """Test that subdividing non-existent edge raises ValueError."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})]
        )
        
        with pytest.raises(ValueError, match="Edge \\(3, 2, key=0\\) not found"):
            _subdivide_edge(net, 3, 2, 0)

    def test_subdivide_edge_preserves_other_edges(self) -> None:
        """Test that subdividing one edge doesn't affect other edges."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 3, 'v': 1, 'key': 0, 'branch_length': 1.0},
                {'u': 3, 'v': 2, 'key': 0, 'branch_length': 2.0},
                (3, 4)  # Node 3 needs degree >= 3
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        graph, subdiv_node = _subdivide_edge(net, 3, 1, 0)
        
        # Other edge should still exist
        assert graph._undirected.has_edge(3, 2, key=0)
        edge_data = graph._undirected[3][2][0]
        assert edge_data.get('branch_length') == 2.0

    def test_subdivide_edge_preserves_nodes(self) -> None:
        """Test that all original nodes are preserved."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'root'})]
        )
        
        graph, subdiv_node = _subdivide_edge(net, 3, 1, 0)
        
        # All original nodes should be present
        assert 1 in graph.nodes()
        assert 3 in graph.nodes()
        assert subdiv_node in graph.nodes()
        
        # Node attributes should be preserved
        assert graph._undirected.nodes[1].get('label') == 'A'
        assert graph._undirected.nodes[3].get('label') == 'root'

