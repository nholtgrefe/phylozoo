"""
Tests for SemiDirectedPhyNetwork and MixedPhyNetwork isomorphism functions.
"""

import pytest

from phylozoo.core.network.sdnetwork import (
    SemiDirectedPhyNetwork,
    MixedPhyNetwork,
    isomorphism,
)
from tests.fixtures.sd_networks import (
    SDTREE_SMALL_BINARY,
    SDTREE_SIMPLE,
    LEVEL_1_SDNETWORK_SINGLE_HYBRID,
)


class TestIsIsomorphic:
    """Tests for is_isomorphic function."""
    
    def test_simple_isomorphic_same_labels(self) -> None:
        """Test that networks with same structure and labels are isomorphic."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2) is True
    
    def test_simple_non_isomorphic_different_labels(self) -> None:
        """Test that networks with different labels are not isomorphic."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'D'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2) is False
    
    def test_same_network(self) -> None:
        """Test that a network is isomorphic to itself."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net, net) is True
    
    def test_fixture_networks_same(self) -> None:
        """Test that fixture networks are isomorphic to themselves."""
        assert isomorphism.is_isomorphic(SDTREE_SMALL_BINARY, SDTREE_SMALL_BINARY) is True
        assert isomorphism.is_isomorphic(SDTREE_SIMPLE, SDTREE_SIMPLE) is True
        assert isomorphism.is_isomorphic(
            LEVEL_1_SDNETWORK_SINGLE_HYBRID,
            LEVEL_1_SDNETWORK_SINGLE_HYBRID
        ) is True
    
    def test_fixture_networks_different(self) -> None:
        """Test that different fixture networks are not isomorphic."""
        # SDTREE_SMALL_BINARY and SDTREE_SIMPLE are actually isomorphic (both star trees with same labels)
        # So test with a tree vs a hybrid network instead
        assert isomorphism.is_isomorphic(
            SDTREE_SMALL_BINARY,
            LEVEL_1_SDNETWORK_SINGLE_HYBRID
        ) is False
    
    def test_with_additional_node_attributes(self) -> None:
        """Test isomorphism with additional node attributes."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A', 'type': 'leaf'}),
                (2, {'label': 'B', 'type': 'leaf'}),
                (4, {'label': 'C', 'type': 'leaf'})
            ]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[
                (6, {'label': 'A', 'type': 'leaf'}),
                (7, {'label': 'B', 'type': 'leaf'}),
                (8, {'label': 'C', 'type': 'leaf'})
            ]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, node_attrs=['type']) is True
    
    def test_with_different_node_attributes(self) -> None:
        """Test that networks with different node attributes are not isomorphic."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A', 'type': 'leaf'}),
                (2, {'label': 'B', 'type': 'leaf'}),
                (4, {'label': 'C', 'type': 'leaf'})
            ]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[
                (6, {'label': 'A', 'type': 'internal'}),
                (7, {'label': 'B', 'type': 'leaf'}),
                (8, {'label': 'C', 'type': 'leaf'})
            ]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, node_attrs=['type']) is False
    
    def test_with_edge_attributes(self) -> None:
        """Test isomorphism with edge attributes."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3},
                (3, 4)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 5, 'v': 6, 'branch_length': 0.5},
                {'u': 5, 'v': 7, 'branch_length': 0.3},
                (5, 8)
            ],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, edge_attrs=['branch_length']) is True
    
    def test_with_different_edge_attributes(self) -> None:
        """Test that networks with different edge attributes are not isomorphic."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3},
                (3, 4)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 5, 'v': 6, 'branch_length': 0.5},
                {'u': 5, 'v': 7, 'branch_length': 0.4},  # Different
                (5, 8)
            ],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(net1, net2, edge_attrs=['branch_length']) is False
    
    def test_with_graph_attributes(self) -> None:
        """Test isomorphism with graph attributes."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            attributes={'version': '1.0', 'source': 'test'}
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})],
            attributes={'version': '1.0', 'source': 'test'}
        )
        
        assert isomorphism.is_isomorphic(
            net1, net2, graph_attrs=['version', 'source']
        ) is True
    
    def test_with_different_graph_attributes(self) -> None:
        """Test that networks with different graph attributes are not isomorphic."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            attributes={'version': '1.0'}
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})],
            attributes={'version': '2.0'}
        )
        
        assert isomorphism.is_isomorphic(net1, net2, graph_attrs=['version']) is False
    
    def test_mixed_network_type(self) -> None:
        """Test that SemiDirectedPhyNetwork and MixedPhyNetwork can be compared."""
        sd_net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        mixed_net = MixedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})]
        )
        
        assert isomorphism.is_isomorphic(sd_net, mixed_net) is True
    
    def test_hybrid_networks(self) -> None:
        """Test isomorphism with hybrid networks."""
        net1 = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)],
            nodes=[
                (3, {'label': 'C'}),
                (7, {'label': 'D'}),
                (1, {'label': 'A'}),
                (2, {'label': 'B'})
            ]
        )
        net2 = SemiDirectedPhyNetwork(
            directed_edges=[(10, 9), (11, 9)],
            undirected_edges=[(10, 12), (10, 11), (11, 13), (9, 14), (14, 15), (14, 16)],
            nodes=[
                (12, {'label': 'C'}),
                (13, {'label': 'D'}),
                (15, {'label': 'A'}),
                (16, {'label': 'B'})
            ]
        )
        
        assert isomorphism.is_isomorphic(net1, net2) is True
    
    def test_label_always_checked(self) -> None:
        """Test that labels are always checked even if node_attrs is None."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'D'})]
        )
        
        # Should be False even without specifying node_attrs
        assert isomorphism.is_isomorphic(net1, net2) is False
    
    def test_label_in_additional_attrs(self) -> None:
        """Test that 'label' is included even when specified in node_attrs."""
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A', 'type': 'leaf'}),
                (2, {'label': 'B', 'type': 'leaf'}),
                (4, {'label': 'C', 'type': 'leaf'})
            ]
        )
        net2 = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 6), (5, 7), (5, 8)],
            nodes=[
                (6, {'label': 'A', 'type': 'leaf'}),
                (7, {'label': 'B', 'type': 'leaf'}),
                (8, {'label': 'C', 'type': 'leaf'})
            ]
        )
        
        # Should work even if 'label' is in node_attrs (no duplicate)
        assert isomorphism.is_isomorphic(net1, net2, node_attrs=['label', 'type']) is True

