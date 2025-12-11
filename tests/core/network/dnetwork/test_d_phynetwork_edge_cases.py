"""
Comprehensive tests for DirectedPhyNetwork edge cases and error scenarios.

This module tests edge cases including:
- Empty networks
- Single node/edge networks
- Very large networks
- High degree nodes
- Many parallel edges
- Deep/wide trees
- Complex topologies
- Invalid structures
- Type errors
"""

import math
import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestEmptyNetworkEdgeCases:
    """Test edge cases for empty networks."""

    def test_empty_network_validation(self) -> None:
        """Test validation of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        with pytest.warns(UserWarning, match="Empty network.*no nodes.*detected"):
            net.validate()

    def test_empty_network_properties(self) -> None:
        """Test properties of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        assert len(net.leaves) == 0
        assert len(net.taxa) == 0
        # Empty network has no root - accessing root_node should raise
        with pytest.raises(ValueError, match="Network has no root node"):
            _ = net.root_node
        # internal_nodes now handles empty networks gracefully
        assert len(net.internal_nodes) == 0

    def test_empty_network_operations(self) -> None:
        """Test operations on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert len(list(net)) == 0
        assert len(net) == 0
        assert 999 not in net


class TestSingleEdgeNetwork:
    """Test edge cases for single edge network."""

    def test_single_edge_network(self) -> None:
        """Test network with single edge."""
        net = DirectedPhyNetwork(edges=[(1, 2)], nodes=[(2, {'label': 'A'})])
        assert net.number_of_nodes() == 2
        assert net.number_of_edges() == 1
        assert net.root_node == 1
        assert net.leaves == {2}


class TestLargeNetworks:
    """Test edge cases with large networks."""

    def test_network_100_nodes(self) -> None:
        """Test network with 100 nodes."""
        # Create star tree with 100 leaves
        edges = [(100, i) for i in range(1, 100)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 100)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.number_of_nodes() == 100
        assert net.number_of_edges() == 99
        assert len(net.leaves) == 99
        net.validate()

    def test_network_100_edges(self) -> None:
        """Test network with 100 edges."""
        # Create a simpler structure: root with many children, each leading to a subtree
        edges = []
        nodes = []
        root = 10000
        
        # Root with 50 children, each child has 2 leaves (50 * 2 = 100 edges)
        for i in range(50):
            child = 20000 + i
            leaf1 = 30000 + 2 * i
            leaf2 = 30000 + 2 * i + 1
            edges.append((root, child))
            edges.append((child, leaf1))
            edges.append((child, leaf2))
            nodes.append((leaf1, {"label": f"Taxon{2*i}"}))
            nodes.append((leaf2, {"label": f"Taxon{2*i+1}"}))
        
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.number_of_edges() == 150  # 50 + 50 + 50 = 150 edges
        net.validate()

    def test_network_10_hybrid_nodes(self) -> None:
        """Test network with 10 hybrid nodes."""
        edges = []
        nodes = []
        root = 10000
        
        # Create 10 independent hybrid events
        for i in range(10):
            tree1 = 1000 + 2 * i
            tree2 = 1000 + 2 * i + 1
            hybrid = 2000 + i
            leaf = 3000 + i
            
            edges.append((root, tree1))
            edges.append((root, tree2))
            edges.append((tree1, hybrid))
            edges.append((tree1, 4000 + i),)  # Tree node 1 also has another child
            edges.append((tree2, hybrid))
            edges.append((tree2, 5000 + i),)  # Tree node 2 also has another child
            edges.append((hybrid, leaf))
            nodes.append((leaf, {"label": f"Taxon{i}"}))
            nodes.append((4000 + i, {"label": f"TaxonA{i}"}))
            nodes.append((5000 + i, {"label": f"TaxonB{i}"}))
        
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert len(net.hybrid_nodes) == 10
        net.validate()


class TestHighDegreeNodes:
    """Test edge cases with high degree nodes."""

    def test_high_outdegree_node(self) -> None:
        """Test node with high out-degree (many children)."""
        # Root with 50 children
        edges = [(1, i) for i in range(2, 52)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(2, 52)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.outdegree(1) == 50
        net.validate()

    def test_high_indegree_hybrid(self) -> None:
        """Test hybrid node with high in-degree (many parents)."""
        # Hybrid with 20 parents
        edges = []
        nodes = []
        root = 10000
        
        for i in range(20):
            tree_node = 20000 + i
            leaf = 30000 + i
            edges.append((root, tree_node))
            edges.append((tree_node, 100))  # All point to hybrid
            edges.append((tree_node, leaf))  # Tree node also has another child
            nodes.append((leaf, {"label": f"Taxon{i}"}))
        
        edges.append((100, 200))  # Hybrid to leaf
        nodes.append((200, {"label": "Leaf"}))
        
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.indegree(100) == 20
        assert 100 in net.hybrid_nodes
        net.validate()


class TestParallelEdges:
    """Test edge cases with parallel edges."""

    def test_many_parallel_edges(self) -> None:
        """Test many parallel edges between same nodes."""
        # 50 parallel edges from 5 to 4 (hybrid node)
        edges = []
        root = 10000
        edges.append((root, 5))
        edges.append((root, 6))
        edges.extend([(5, 4, i) for i in range(50)])  # 50 parallel edges
        edges.append((5, 8))  # Tree node 5 also has another child
        edges.append((6, 4))
        edges.append((6, 9))  # Tree node 6 also has another child
        edges.append((4, 2))  # Hybrid to leaf
        net = DirectedPhyNetwork(edges=edges, nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})])
        # Total edges: 2 (root->5, root->6) + 50 (parallel) + 1 (5->8) + 1 (6->4) + 1 (6->9) + 1 (4->2) = 56
        assert net.number_of_edges() == 56
        assert net.indegree(4) == 51  # 50 from 5, 1 from 6
        net.validate()

    def test_parallel_edges_to_hybrid_with_gamma(self) -> None:
        """Test parallel edges to hybrid with gamma values."""
        edges = [
            (7, 5), (7, 6),  # Root to tree nodes
            (5, 4, 0), (5, 4, 1),  # Parallel edges (no gamma - valid)
            (5, 8),  # Tree node 5 also has another child
            (6, 4),
            (6, 9),  # Tree node 6 also has another child
            (4, 2)  # Hybrid to leaf
        ]
        net = DirectedPhyNetwork(edges=edges, nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})])
        # No gamma values set, so validation should pass
        net.validate()


class TestDeepTrees:
    """Test edge cases with very deep trees."""

    def test_very_deep_tree(self) -> None:
        """Test very deep tree (chain)."""
        edges = []
        nodes = []
        node_id = 10000  # Start with high ID to avoid conflicts
        
        def build_deep_tree(parent, level, max_level):
            nonlocal node_id
            if level >= max_level:
                # At max level, add 2 leaves to ensure parent has out-degree >= 2
                leaf1 = node_id
                node_id += 1
                leaf2 = node_id
                node_id += 1
                edges.append((parent, leaf1))
                edges.append((parent, leaf2))
                nodes.append((leaf1, {"label": f"Taxon{leaf1}"}))
                nodes.append((leaf2, {"label": f"Taxon{leaf2}"}))
                return
            
            left = node_id
            node_id += 1
            right = node_id
            node_id += 1
            
            edges.append((parent, left))
            edges.append((parent, right))
            build_deep_tree(left, level + 1, max_level)
            build_deep_tree(right, level + 1, max_level)
        
        root = 1
        build_deep_tree(root, 0, 6)  # 6 levels gives deep tree
        
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.number_of_nodes() >= 50
        assert net.number_of_edges() >= 50
        net.validate()
        assert net.root_node == root

    def test_deep_binary_tree(self) -> None:
        """Test deep binary tree."""
        # Binary tree with 10 levels (2^10 = 1024 nodes max)
        edges = []
        nodes = []
        node_id = 10000
        
        def build_tree(parent, level, max_level):
            nonlocal node_id
            if level >= max_level:
                # At max level, add 2 leaves to ensure parent has out-degree >= 2
                leaf1 = node_id
                node_id += 1
                leaf2 = node_id
                node_id += 1
                edges.append((parent, leaf1))
                edges.append((parent, leaf2))
                nodes.append((leaf1, {"label": f"Taxon{leaf1}"}))
                nodes.append((leaf2, {"label": f"Taxon{leaf2}"}))
                return
            
            left = node_id
            node_id += 1
            right = node_id
            node_id += 1
            
            edges.append((parent, left))
            edges.append((parent, right))
            build_tree(left, level + 1, max_level)
            build_tree(right, level + 1, max_level)
        
        root = 1
        build_tree(root, 0, 6)  # 6 levels
        
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        net.validate()
        assert net.root_node == root


class TestWideTrees:
    """Test edge cases with very wide trees."""

    def test_very_wide_tree(self) -> None:
        """Test very wide tree (many leaves from root)."""
        # Root with 200 leaves
        edges = [(1, i) for i in range(2, 202)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(2, 202)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.number_of_nodes() == 201
        assert len(net.leaves) == 200
        net.validate()


class TestComplexTopologies:
    """Test edge cases with complex network topologies."""

    def test_multiple_levels_hybridization(self) -> None:
        """Test network with multiple levels of hybridization."""
        # Level-1: hybrid 4
        # Level-2: hybrid 5 (parent of hybrid 4)
        # Need to ensure all tree nodes have out-degree >= 2
        edges = [
            (10, 7), (10, 8),  # Root splits
            (7, 5), (7, 6),    # Tree node 7 splits (5 and 6 are children)
            (7, 11),           # Tree node 7 also has another child
            (8, 5), (8, 9),    # Tree node 8 splits (5 is hybrid, 9 is tree node)
            (8, 13),           # Tree node 8 also has another child
            (5, 4), (6, 4),    # Hybrid 4 (5 and 6 point to it)
            (6, 12),           # Tree node 6 also has another child
            (9, 2), (9, 14),   # Tree node 9 splits (out-degree 2)
            (4, 1)             # Hybrid to leaf
        ]
        net = DirectedPhyNetwork(edges=edges, nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'}), (13, {'label': 'E'}), (14, {'label': 'F'})])
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes
        net.validate()

    def test_nested_hybridization(self) -> None:
        """Test nested hybridization."""
        # Hybrid 4 is child of hybrid 5
        # Need to ensure all tree nodes have out-degree >= 2
        edges = [
            (10, 7), (10, 8),
            (7, 5), (7, 6), (7, 11),  # Tree node 7 splits
            (8, 5), (8, 9), (8, 12),  # Tree node 8 splits
            (5, 4), (6, 4),           # Hybrid 4 (5 and 6 point to it)
            (6, 13),                  # Tree node 6 also has another child
            (9, 2), (9, 14),          # Tree node 9 splits (out-degree 2)
            (4, 1)                    # Hybrid to leaf
        ]
        net = DirectedPhyNetwork(edges=edges, nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'}), (13, {'label': 'E'}), (14, {'label': 'F'})])
        # Both 4 and 5 are hybrids
        assert len(net.hybrid_nodes) == 2
        net.validate()

    def test_multiple_independent_hybrids(self) -> None:
        """Test multiple independent hybrid events."""
        edges = []
        nodes = []
        root = 10000
        
        # 5 independent hybrid events
        for i in range(5):
            t1 = 1000 + 2 * i
            t2 = 1000 + 2 * i + 1
            h = 2000 + i
            l = 3000 + i
            
            edges.append((root, t1))
            edges.append((root, t2))
            edges.append((t1, h))
            edges.append((t1, 4000 + i))  # Tree node t1 also has another child
            edges.append((t2, h))
            edges.append((t2, 5000 + i))  # Tree node t2 also has another child
            edges.append((h, l))
            nodes.append((l, {"label": f"Taxon{i}"}))
            nodes.append((4000 + i, {"label": f"TaxonA{i}"}))
            nodes.append((5000 + i, {"label": f"TaxonB{i}"}))
        
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert len(net.hybrid_nodes) == 5
        net.validate()


class TestInvalidStructures:
    """Test invalid network structures (should raise errors)."""

    def test_cycle_simple(self) -> None:
        """Test simple cycle detection."""
        # Cycle: 1 -> 2 -> 3 -> 1
        # Need to ensure nodes in cycle are not leaves (they have outgoing edges)
        with pytest.raises(ValueError, match="directed cycles"):
            DirectedPhyNetwork(
                edges=[(1, 2), (2, 3), (3, 1), (1, 4)],  # Add leaf 4 from node 1
                nodes=[(4, {'label': 'A'})]  # Only leaf 4 is in taxa
            )

    def test_cycle_complex(self) -> None:
        """Test complex cycle detection."""
        # Cycle: 2->3->4->2
        # Need to ensure nodes in cycle are not leaves
        with pytest.raises(ValueError, match="directed cycles"):
            DirectedPhyNetwork(
                edges=[(1, 2), (2, 3), (3, 4), (4, 2), (1, 5)],  # Add leaf 5 from node 1
                nodes=[(5, {'label': 'A'})]  # Only leaf 5 is in taxa
            )

    def test_multiple_roots(self) -> None:
        """Test multiple roots detection."""
        with pytest.raises(ValueError, match="multiple root nodes"):
            DirectedPhyNetwork(
                edges=[(1, 3), (2, 3)],
                nodes=[(3, {'label': 'A'})]
            )

    def test_leaf_wrong_indegree(self) -> None:
        """Test leaf with wrong in-degree."""
        # Leaf 3 has in-degree 2 (from 1 and 2)
        # Need single root, so add root that points to 1 and 2
        with pytest.raises(ValueError, match="in-degree"):
            DirectedPhyNetwork(
                edges=[(10, 1), (10, 2), (1, 3), (2, 3)],  # Root 10, nodes 1 and 2 both point to 3
                nodes=[(3, {'label': 'A'})]
            )

    def test_internal_node_invalid_degrees(self) -> None:
        """Test internal node with invalid degrees."""
        # Node 2 has in-degree 1, out-degree 1 (invalid)
        with pytest.raises(ValueError, match="Internal node"):
            DirectedPhyNetwork(
                edges=[(1, 2), (2, 3)],
                nodes=[(3, {'label': 'A'})]
            )


class TestInvalidAttributes:
    """Test invalid attribute values (should raise errors)."""

    def test_bootstrap_out_of_range_negative(self) -> None:
        """Test bootstrap < 0.0."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': -0.1}],
                nodes=[(1, {'label': 'A'})]
            )

    def test_bootstrap_out_of_range_positive(self) -> None:
        """Test bootstrap > 1.0."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': 1.1}],
                nodes=[(1, {'label': 'A'})]
            )

    def test_gamma_out_of_range_negative(self) -> None:
        """Test gamma < 0.0."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': -0.1},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 1.1},
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
            )

    def test_gamma_sum_not_one(self) -> None:
        """Test gamma sum != 1.0."""
        with pytest.raises(ValueError, match="must sum to exactly 1.0"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 0.3},  # Sum = 0.9
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
            )

    def test_gamma_partial_specification(self) -> None:
        """Test partial gamma specification (should fail)."""
        with pytest.raises(ValueError, match="ALL incoming edges must have gamma values"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4},  # Missing gamma
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
            )


class TestTypeErrors:
    """Test type errors in attributes."""

    def test_bootstrap_non_numeric(self) -> None:
        """Test non-numeric bootstrap."""
        with pytest.raises(ValueError, match="must be numeric"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': 'invalid'}],
                nodes=[(1, {'label': 'A'})]
            )

    def test_gamma_non_numeric(self) -> None:
        """Test non-numeric gamma."""
        with pytest.raises(ValueError, match="must be numeric"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 'invalid'},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 1.0},
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
            )

    def test_bootstrap_list_type(self) -> None:
        """Test bootstrap as list (wrong type)."""
        with pytest.raises(ValueError, match="must be numeric"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': [0.5]}],
                nodes=[(1, {'label': 'A'})]
            )


class TestBoundaryValues:
    """Test boundary values for attributes."""

    def test_bootstrap_exactly_zero(self) -> None:
        """Test bootstrap exactly 0.0."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'bootstrap': 0.0}],
            nodes=[(1, {'label': 'A'})]
        )
        assert net.get_bootstrap(3, 1) == 0.0
        net.validate()

    def test_bootstrap_exactly_one(self) -> None:
        """Test bootstrap exactly 1.0."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'bootstrap': 1.0}],
            nodes=[(1, {'label': 'A'})]
        )
        assert net.get_bootstrap(3, 1) == 1.0
        net.validate()

    def test_gamma_exactly_zero(self) -> None:
        """Test gamma exactly 0.0."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.0},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'gamma': 1.0},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.get_gamma(5, 4) == 0.0
        net.validate()

    def test_gamma_exactly_one(self) -> None:
        """Test gamma exactly 1.0."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 1.0},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'gamma': 0.0},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.get_gamma(5, 4) == 1.0
        net.validate()

    def test_gamma_sum_floating_point_precision(self) -> None:
        """Test gamma sum with floating point precision issues."""
        # Values that sum to 1.0 but might have precision issues
        net = DirectedPhyNetwork(
            edges=[
                (20, 5), (20, 6), (20, 7),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 1.0 / 3.0},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'gamma': 1.0 / 3.0},
                {'u': 6, 'v': 9},
                {'u': 7, 'v': 4, 'gamma': 1.0 / 3.0},
                {'u': 7, 'v': 10},
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})]
        )
        # Should pass with tolerance
        net.validate()


class TestSpecialValues:
    """Test special values (NaN, infinity) in attributes."""

    def test_nan_branch_length(self) -> None:
        """Test NaN branch length (should be allowed, not validated)."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': float('nan')}],
            nodes=[(1, {'label': 'A'})]
        )
        result = net.get_branch_length(3, 1)
        assert result is not None
        assert math.isnan(result)

    def test_infinity_branch_length(self) -> None:
        """Test infinity branch length (should be allowed, not validated)."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': float('inf')}],
            nodes=[(1, {'label': 'A'})]
        )
        result = net.get_branch_length(3, 1)
        assert result is not None
        assert math.isinf(result)

    def test_nan_bootstrap(self) -> None:
        """Test NaN bootstrap (should fail validation)."""
        # NaN is not in [0.0, 1.0], so should fail
        # The validation now explicitly checks for NaN using math.isnan()
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': float('nan')}],
                nodes=[(1, {'label': 'A'})]
            )

    def test_infinity_bootstrap(self) -> None:
        """Test infinity bootstrap (should fail validation)."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': float('inf')}],
                nodes=[(1, {'label': 'A'})]
            )

