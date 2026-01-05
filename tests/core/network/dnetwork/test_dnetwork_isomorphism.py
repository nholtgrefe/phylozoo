"""
Tests for DirectedPhyNetwork isomorphism functions.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork, isomorphism
from tests.fixtures.directed_networks import (
    DTREE_SMALL_BINARY,
    DTREE_SIMPLE,
    LEVEL_1_DNETWORK_SINGLE_HYBRID,
)


class TestIsIsomorphic:
    """Tests for is_isomorphic function."""
    
    def test_simple_isomorphic_same_labels(self) -> None:
        """Test that networks with same structure and labels are isomorphic."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2) is True
    
    def test_simple_non_isomorphic_different_labels(self) -> None:
        """Test that networks with different labels are not isomorphic."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'D'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2) is False
    
    def test_same_network(self) -> None:
        """Test that a network is isomorphic to itself."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net, net) is True
    
    def test_fixture_networks_same(self) -> None:
        """Test that fixture networks are isomorphic to themselves."""
        assert isomorphism.is_isomorphic(DTREE_SMALL_BINARY, DTREE_SMALL_BINARY) is True
        assert isomorphism.is_isomorphic(DTREE_SIMPLE, DTREE_SIMPLE) is True
        assert isomorphism.is_isomorphic(
            LEVEL_1_DNETWORK_SINGLE_HYBRID,
            LEVEL_1_DNETWORK_SINGLE_HYBRID
        ) is True
    
    def test_fixture_networks_different(self) -> None:
        """Test that different fixture networks are not isomorphic."""
        assert isomorphism.is_isomorphic(DTREE_SMALL_BINARY, DTREE_SIMPLE) is False
        assert isomorphism.is_isomorphic(
            DTREE_SMALL_BINARY,
            LEVEL_1_DNETWORK_SINGLE_HYBRID
        ) is False
    
    def test_with_additional_node_attributes(self) -> None:
        """Test isomorphism with additional node attributes."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A', 'type': 'leaf'}),
                (2, {'label': 'B', 'type': 'leaf'}),
                (4, {'label': 'C', 'type': 'leaf'})
            ]
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[
                (6, {'label': 'A', 'type': 'leaf'}),
                (7, {'label': 'B', 'type': 'leaf'}),
                (8, {'label': 'C', 'type': 'leaf'})
            ]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, node_attrs=['type']) is True
    
    def test_with_different_node_attributes(self) -> None:
        """Test that networks with different node attributes are not isomorphic."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A', 'type': 'leaf'}),
                (2, {'label': 'B', 'type': 'leaf'}),
                (4, {'label': 'C', 'type': 'leaf'})
            ]
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[
                (6, {'label': 'A', 'type': 'internal'}),
                (7, {'label': 'B', 'type': 'leaf'}),
                (8, {'label': 'C', 'type': 'leaf'})
            ]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, node_attrs=['type']) is False
    
    def test_with_edge_attributes(self) -> None:
        """Test isomorphism with edge attributes."""
        net1 = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3},
                (3, 4)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 6, 'branch_length': 0.5},
                {'u': 5, 'v': 7, 'branch_length': 0.3},
                (5, 8)
            ],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, edge_attrs=['branch_length']) is True
    
    def test_with_different_edge_attributes(self) -> None:
        """Test that networks with different edge attributes are not isomorphic."""
        net1 = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3},
                (3, 4)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 6, 'branch_length': 0.5},
                {'u': 5, 'v': 7, 'branch_length': 0.4},  # Different
                (5, 8)
            ],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, edge_attrs=['branch_length']) is False
    
    def test_with_graph_attributes(self) -> None:
        """Test isomorphism with graph attributes."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            attributes={'version': '1.0', 'source': 'test'}
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})],
            attributes={'version': '1.0', 'source': 'test'}
        )
        
        assert isomorphism.is_isomorphic(
            net1, net2, graph_attrs=['version', 'source']
        ) is True
    
    def test_with_different_graph_attributes(self) -> None:
        """Test that networks with different graph attributes are not isomorphic."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            attributes={'version': '1.0'}
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})],
            attributes={'version': '2.0'}
        )
        
        assert isomorphism.is_isomorphic(net1, net2, graph_attrs=['version']) is False
    
    def test_hybrid_networks(self) -> None:
        """Test isomorphism with hybrid networks."""
        net1 = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4), (6, 4),
                (5, 8), (6, 9),
                (4, 1)
            ],
            nodes=[
                (1, {'label': 'A'}),
                (8, {'label': 'B'}),
                (9, {'label': 'C'})
            ]
        )
        net2 = DirectedPhyNetwork(
            edges=[
                (12, 10), (12, 11),
                (10, 9), (11, 9),
                (10, 13), (11, 14),
                (9, 15)
            ],
            nodes=[
                (15, {'label': 'A'}),
                (13, {'label': 'B'}),
                (14, {'label': 'C'})
            ]
        )
        
        assert isomorphism.is_isomorphic(net1, net2) is True
    
    def test_label_always_checked(self) -> None:
        """Test that labels are always checked even if node_attrs is None."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'D'})]
        )
        
        # Should be False even without specifying node_attrs
        assert isomorphism.is_isomorphic(net1, net2) is False
    
    def test_label_in_additional_attrs(self) -> None:
        """Test that 'label' is included even when specified in node_attrs."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A', 'type': 'leaf'}),
                (2, {'label': 'B', 'type': 'leaf'}),
                (4, {'label': 'C', 'type': 'leaf'})
            ]
        )
        net2 = DirectedPhyNetwork(
            edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[
                (6, {'label': 'A', 'type': 'leaf'}),
                (7, {'label': 'B', 'type': 'leaf'}),
                (8, {'label': 'C', 'type': 'leaf'})
            ]
        )
        
        # Should work even if 'label' is in node_attrs (no duplicate)
        assert isomorphism.is_isomorphic(net1, net2, node_attrs=['label', 'type']) is True

