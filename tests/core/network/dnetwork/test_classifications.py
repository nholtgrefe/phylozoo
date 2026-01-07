"""
Tests for classification functions for directed phylogenetic networks.

This module tests classification functions including:
- level
- vertex_level
- reticulation_number
- is_binary
- is_tree
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import (
    is_binary,
    is_simple,
    is_stackfree,
    is_tree,
    level,
    reticulation_number,
    vertex_level,
)


class TestLevel:
    """Test cases for level() function."""

    def test_level_empty_network(self) -> None:
        """Test level in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert level(net) == 0

    def test_level_tree(self) -> None:
        """Test level in tree (no hybrids)."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert level(net) == 0

    def test_level_single_hybrid(self) -> None:
        """Test level with single hybrid node."""
        # Network: root -> tree nodes -> hybrid -> leaf
        # Hybrid node 4 has 2 incoming edges (from 5 and 6)
        # Blob contains {4, 5, 6, 8}
        # Hybrid edges in blob: (5, 4), (6, 4) = 2
        # Hybrid nodes in blob: {4} = 1
        # Level = 2 - 1 = 1
        net = DirectedPhyNetwork(
            edges=[(8, 5), (8, 6), (5, 4), (5, 1), (6, 4), (6, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        assert level(net) == 1

    def test_level_multiple_blobs(self) -> None:
        """Test level with multiple blobs."""
        # Network with two hybrid nodes in the same blob (bi-edge connected via node 10)
        # With bi-edge connectivity, both hybrid regions are in one blob since there are no bridges
        # Single blob contains both hybrid nodes 4 and 7
        # Hybrid edges: (5, 4), (6, 4), (8, 7), (9, 7) = 4 edges
        # Hybrid nodes: 4, 7 = 2 nodes
        # Level = 4 - 2 = 2
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6), (10, 8), (10, 9),
                (5, 4), (6, 4), (5, 1), (6, 2),
                (8, 7), (9, 7), (8, 3), (9, 11),
                (4, 12), (7, 13)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'})
            ]
        )
        assert level(net) == 2

    def test_level_nested_hybrids(self) -> None:
        """Test level with nested hybrids in same blob."""
        # Network with nested hybrids in same blob
        # Hybrid 5: edges (7, 5), (8, 5) = 2 edges, 1 node
        # Hybrid 4: edges (5, 4), (6, 4) = 2 edges, 1 node
        # Total in blob: 4 edges, 2 nodes -> level = 4 - 2 = 2
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8), (10, 6),
                (7, 5), (8, 5), (7, 1), (8, 2),
                (5, 4), (6, 4), (6, 3),
                (4, 9)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (9, {'label': 'D'})
            ]
        )
        assert level(net) == 2


class TestVertexLevel:
    """Test cases for vertex_level() function."""

    def test_vertex_level_empty_network(self) -> None:
        """Test vertex_level in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert vertex_level(net) == 0

    def test_vertex_level_tree(self) -> None:
        """Test vertex_level in tree (no hybrids)."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert vertex_level(net) == 0

    def test_vertex_level_single_hybrid(self) -> None:
        """Test vertex_level with single hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(8, 5), (8, 6), (5, 4), (5, 1), (6, 4), (6, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        # Blob contains hybrid node 4
        assert vertex_level(net) == 1

    def test_vertex_level_multiple_hybrids_same_blob(self) -> None:
        """Test vertex_level with multiple hybrids in same blob."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8), (10, 6),
                (7, 5), (8, 5), (7, 1), (8, 2),
                (5, 4), (6, 4), (6, 3),
                (4, 9)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (9, {'label': 'D'})
            ]
        )
        # Blob contains hybrid nodes 4 and 5
        assert vertex_level(net) == 2

    def test_vertex_level_multiple_blobs(self) -> None:
        """Test vertex_level with multiple blobs."""
        # Network with two hybrid nodes in the same blob (bi-edge connected via node 10)
        # With bi-edge connectivity, both hybrid regions are in one blob since there are no bridges
        # Single blob contains both hybrid nodes 4 and 7
        # Vertex level = max hybrid nodes in blob = 2
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6), (10, 8), (10, 9),
                (5, 4), (6, 4), (5, 1), (6, 2),
                (8, 7), (9, 7), (8, 3), (9, 11),
                (4, 12), (7, 13)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'})
            ]
        )
        # Single blob has 2 hybrid nodes, maximum is 2
        assert vertex_level(net) == 2


class TestReticulationNumber:
    """Test cases for reticulation_number() function."""

    def test_reticulation_number_empty_network(self) -> None:
        """Test reticulation_number in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert reticulation_number(net) == 0

    def test_reticulation_number_tree(self) -> None:
        """Test reticulation_number in tree (no hybrids)."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert reticulation_number(net) == 0

    def test_reticulation_number_single_hybrid(self) -> None:
        """Test reticulation_number with single hybrid node."""
        # Hybrid node 4 has 2 incoming edges
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # 2 hybrid edges - 1 hybrid node = 1
        assert reticulation_number(net) == 1

    def test_reticulation_number_multiple_hybrids(self) -> None:
        """Test reticulation_number with multiple hybrid nodes."""
        # Hybrid node 4: 2 edges, Hybrid node 5: 2 edges
        # Node 6 needs out-degree >= 2, so add another child
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8),
                (7, 5), (7, 6),
                (8, 5), (8, 9),
                (5, 4), (6, 4), (6, 11),
                (4, 1), (9, 2), (9, 3)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (11, {'label': 'D'})
            ]
        )
        # 4 hybrid edges - 2 hybrid nodes = 2
        assert reticulation_number(net) == 2


class TestIsBinary:
    """Test cases for is_binary() function."""

    def test_is_binary_empty_network(self) -> None:
        """Test is_binary in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert is_binary(net) is True

    def test_is_binary_simple_tree(self) -> None:
        """Test is_binary in simple binary tree."""
        # Root has degree 2, internal node has degree 3
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'})]
        )
        assert is_binary(net) is True

    def test_is_binary_non_binary_root(self) -> None:
        """Test is_binary with root having wrong degree."""
        # Root has degree 3, should be 2
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert is_binary(net) is False

    def test_is_binary_non_binary_internal(self) -> None:
        """Test is_binary with internal node having wrong degree."""
        # Internal node 3 has degree 4, should be 3
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (3, 5), (3, 6), (4, 7)],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}),
                (6, {'label': 'D'}), (7, {'label': 'E'})
            ]
        )
        assert is_binary(net) is False

    def test_is_binary_with_hybrids(self) -> None:
        """Test is_binary in network with hybrid nodes."""
        # Hybrid node 4 has degree 3 (in-degree 2, out-degree 1)
        # Tree nodes 5, 6 have degree 3
        # Root 7 has degree 2
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert is_binary(net) is True

    def test_is_binary_single_node(self) -> None:
        """Test is_binary in single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[], nodes=[(1, {'label': 'A'})])
        assert is_binary(net) is True


class TestIsTree:
    """Test cases for is_tree() function."""

    def test_is_tree_empty_network(self) -> None:
        """Test is_tree in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert is_tree(net) is True

    def test_is_tree_simple_tree(self) -> None:
        """Test is_tree in simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert is_tree(net) is True

    def test_is_tree_with_hybrids(self) -> None:
        """Test is_tree in network with hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert is_tree(net) is False

    def test_is_tree_large_tree(self) -> None:
        """Test is_tree in large binary tree."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 3), (5, 4),
                (6, 1), (6, 2)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'})
            ]
        )
        assert is_tree(net) is True


class TestIsSimple:
    """Test cases for is_simple() function."""

    def test_is_simple_empty_network(self) -> None:
        """Test is_simple in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert is_simple(net) is True

    def test_is_simple_tree_single_internal(self) -> None:
        """Test is_simple in tree with single internal node."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        # Tree with one internal node has one non-leaf blob, so it is simple
        assert is_simple(net) is True

    def test_is_simple_tree_multiple_internal(self) -> None:
        """Test is_simple in tree with multiple internal nodes."""
        net = DirectedPhyNetwork(
            edges=[(5, 3), (5, 4), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        # Tree with multiple internal nodes has multiple non-leaf blobs, so it is not simple
        assert is_simple(net) is False

    def test_is_simple_single_hybrid(self) -> None:
        """Test is_simple with single hybrid node (one non-leaf blob)."""
        # Network: root -> tree nodes -> hybrid -> leaf
        net = DirectedPhyNetwork(
            edges=[(8, 5), (8, 6), (5, 4), (5, 1), (6, 4), (6, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        # Single hybrid creates one non-leaf blob
        assert is_simple(net) is True

    def test_is_simple_multiple_blobs(self) -> None:
        """Test is_simple with multiple non-leaf blobs."""
        # Use a fixture network that has multiple blobs
        from tests.fixtures import directed_networks as dn
        # LEVEL_1_DNETWORK_TWO_BLOBS has two separate hybrid regions
        assert is_simple(dn.LEVEL_1_DNETWORK_TWO_BLOBS) is False

    def test_is_simple_with_fixtures(self) -> None:
        """Test is_simple using example networks from fixtures."""
        from tests.fixtures import directed_networks as dn
        
        # Networks with multiple blobs should not be simple
        assert is_simple(dn.LEVEL_1_DNETWORK_TWO_BLOBS) is False
        assert is_simple(dn.LEVEL_2_DNETWORK_THREE_BLOBS) is False
        assert is_simple(dn.LEVEL_1_DNETWORK_FIVE_BLOBS) is False


class TestIsStackfreeWithFixtures:
    """Test cases for is_stackfree() function using fixture networks."""
    
    def test_is_stackfree_tree_fixtures(self) -> None:
        """Test is_stackfree on tree fixtures (should all be stack-free)."""
        from tests.fixtures import directed_networks as dn
        
        # Trees have no hybrids, so they are stack-free
        assert is_stackfree(dn.DTREE_EMPTY)
        assert is_stackfree(dn.DTREE_SINGLE_NODE)
        assert is_stackfree(dn.DTREE_SMALL_BINARY)
        assert is_stackfree(dn.DTREE_NON_BINARY_SMALL)
    
    def test_is_stackfree_single_hybrid_fixtures(self) -> None:
        """Test is_stackfree on networks with single hybrid (should be stack-free)."""
        from tests.fixtures import directed_networks as dn
        
        # Single hybrid networks have no stacked hybrids
        assert is_stackfree(dn.LEVEL_1_DNETWORK_SINGLE_HYBRID)
        assert is_stackfree(dn.LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY)
    
    def test_is_stackfree_multiple_hybrids_separate_fixtures(self) -> None:
        """Test is_stackfree on networks with multiple hybrids in separate blobs."""
        from tests.fixtures import directed_networks as dn
        
        # Hybrids in separate blobs don't stack
        assert is_stackfree(dn.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE)
    
    def test_is_stackfree_nested_hybrids_fixture(self) -> None:
        """Test is_stackfree on network with nested hybrids (should be stack-free)."""
        from tests.fixtures import directed_networks as dn
        
        # LEVEL_2_DNETWORK_NESTED_HYBRIDS has nested hybrids but they don't stack
        # (each hybrid has a tree node child, not another hybrid)
        assert is_stackfree(dn.LEVEL_2_DNETWORK_NESTED_HYBRIDS)
    
    def test_is_stackfree_chain_hybrids_fixture(self) -> None:
        """Test is_stackfree on network with chain of hybrids (should be stack-free)."""
        from tests.fixtures import directed_networks as dn
        
        # LEVEL_3_DNETWORK_CHAIN_HYBRIDS has a chain but hybrids don't stack
        # (each hybrid has a tree node child)
        assert is_stackfree(dn.LEVEL_3_DNETWORK_CHAIN_HYBRIDS)
    
    def test_is_stackfree_diamond_hybrid_fixture(self) -> None:
        """Test is_stackfree on diamond hybrid structure (should be stack-free)."""
        from tests.fixtures import directed_networks as dn
        
        # LEVEL_2_DNETWORK_DIAMOND_HYBRID has multiple hybrids but they don't stack
        assert is_stackfree(dn.LEVEL_2_DNETWORK_DIAMOND_HYBRID)


class TestIsStackfreeNotStackfree:
    """Test cases for is_stackfree() function on networks with stacked hybrids."""
    
    def test_is_stackfree_stacked_hybrids_simple(self) -> None:
        """Test is_stackfree on network with simple stacked hybrids."""
        # Network: root -> tree nodes -> hybrid 4 -> hybrid 7 -> leaf
        net = DirectedPhyNetwork(
            edges=[
                (9, 5), (9, 6),  # Root to tree nodes
                (5, 4), (5, 8),  # Tree node 5 to hybrid 4 and leaf 8
                (6, 4), (6, 10),  # Tree node 6 to hybrid 4 and leaf 10
                (4, 7), (9, 11), (11, 7), (11, 12),  # Hybrid 4 and tree node 11 lead to hybrid 7
                (7, 1)  # Hybrid 7 to leaf
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (12, {'label': 'D'})]
        )
        assert not is_stackfree(net)
    
    def test_is_stackfree_multiple_stacked_hybrids(self) -> None:
        """Test is_stackfree on network with multiple stacked hybrid pairs."""
        # Network: hybrid 4 -> hybrid 7 -> hybrid 9 -> leaf
        net = DirectedPhyNetwork(
            edges=[
                (11, 5), (11, 6),  # Root to tree nodes
                (5, 4), (5, 12),  # Tree node 5 to hybrid 4 and leaf 12
                (6, 4), (6, 13),  # Tree node 6 to hybrid 4 and leaf 13
                (4, 7), (11, 14), (14, 7), (14, 15),  # Hybrid 4 and tree node 14 lead to hybrid 7
                (7, 9), (11, 16), (16, 9), (16, 17),  # Hybrid 7 and tree node 16 lead to hybrid 9
                (9, 1)  # Hybrid 9 to leaf
            ],
            nodes=[(1, {'label': 'A'}), (12, {'label': 'B'}), (13, {'label': 'C'}), (15, {'label': 'D'}), (17, {'label': 'E'})]
        )
        assert not is_stackfree(net)
    
    def test_is_stackfree_mixed_stacked_and_non_stacked(self) -> None:
        """Test is_stackfree on network with both stacked and non-stacked hybrids."""
        # Network: hybrid 4 -> hybrid 7 (stacked), but hybrid 9 is not stacked
        net = DirectedPhyNetwork(
            edges=[
                (11, 5), (11, 6), (11, 14),  # Root to tree nodes
                (5, 4), (5, 12),  # Tree node 5 to hybrid 4 and leaf 12
                (6, 4), (6, 13),  # Tree node 6 to hybrid 4 and leaf 13
                (4, 7), (14, 7), (14, 15),  # Hybrid 4 and tree node 14 lead to hybrid 7 (stacked)
                (7, 1),  # Hybrid 7 to leaf
                (14, 9), (11, 18), (18, 9), (18, 19),  # Tree node 14 and tree node 18 lead to hybrid 9 (not stacked)
                (9, 2)  # Hybrid 9 to leaf
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (12, {'label': 'C'}), (13, {'label': 'D'}), (15, {'label': 'E'}), (19, {'label': 'F'})]
        )
        assert not is_stackfree(net)  # Has at least one stack

