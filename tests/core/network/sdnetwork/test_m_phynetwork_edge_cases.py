"""
Comprehensive tests for MixedPhyNetwork edge cases and error scenarios.

This module tests edge cases including:
- Empty networks
- Single node/edge networks
- Very large networks
- High degree nodes
- Many parallel edges
- Deep/wide trees
- Complex topologies
- Invalid structures
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestEmptyNetworkEdgeCases:
    """Test edge cases for empty networks."""

    def test_empty_network_validation(self) -> None:
        """Test validation of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert net.validate() is True

    def test_empty_network_properties(self) -> None:
        """Test properties of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        assert len(net.leaves) == 0
        assert len(net.taxa) == 0
        assert len(net.internal_nodes) == 0



class TestSingleEdgeNetwork:
    """Test edge cases for single edge network."""

    def test_single_undirected_edge_network(self) -> None:
        """Test network with single undirected edge."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(1, 2)],
            taxa={2: "A"}
            )
        assert net.number_of_nodes() == 2
        assert net.number_of_edges() == 1
        assert 2 in net.leaves

    def test_single_directed_edge_network(self) -> None:
        """Test network with hybrid node (requires at least 2 incoming directed edges)."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(1, 3), (2, 3)],  # Hybrid node 3
            undirected_edges=[(3, 4), (1, 5), (1, 6), (2, 7), (2, 8)],
            taxa={4: "A", 5: "B", 6: "C", 7: "D", 8: "E"}
            )
        assert net.number_of_nodes() == 8
        assert net.number_of_edges() == 7
        assert 3 in net.hybrid_nodes


class TestLargeNetworks:
    """Test edge cases with large networks."""

    def test_network_100_nodes(self) -> None:
        """Test network with 100 nodes."""
        # Create star tree with 100 leaves
        edges = [(100, i) for i in range(1, 100)]
        taxa = {i: f"Taxon{i}" for i in range(1, 100)}
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.number_of_nodes() == 100
        assert net.number_of_edges() == 99
        assert len(net.leaves) == 99
        with expect_mixed_network_warning():
            assert net.validate() is True

    def test_network_many_hybrids(self) -> None:
        """Test network with many hybrid nodes."""
        # Create network with 10 hybrid nodes, all connected
        directed_edges = []
        undirected_edges = []
        taxa = {}
        
        # Create 10 hybrid events, connected via a central node
        central = 0
        for i in range(10):
            parent1 = 1000 + 2 * i
            parent2 = 1000 + 2 * i + 1
            hybrid = 2000 + i
            leaf = 3000 + i
            
            directed_edges.append((parent1, hybrid))
            directed_edges.append((parent2, hybrid))
            undirected_edges.append((hybrid, leaf))
            undirected_edges.append((parent1, 4000 + 2 * i))
            undirected_edges.append((parent2, 4000 + 2 * i + 1))
            # Connect both parents to central node to ensure degree >= 3
            undirected_edges.append((central, parent1))
            undirected_edges.append((central, parent2))
            taxa[leaf] = f"Taxon{i}"
            taxa[4000 + 2 * i] = f"Taxon{10 + 2 * i}"
            taxa[4000 + 2 * i + 1] = f"Taxon{10 + 2 * i + 1}"
        
        # Central node needs degree >= 3
        undirected_edges.append((central, 5000))
        undirected_edges.append((central, 5001))
        taxa[5000] = "TaxonCentral1"
        taxa[5001] = "TaxonCentral2"
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert len(net.hybrid_nodes) == 10
        with expect_mixed_network_warning():
            assert net.validate() is True


class TestHighDegreeNodes:
    """Test edge cases with high degree nodes."""

    def test_node_degree_10(self) -> None:
        """Test node with degree 10."""
        # Exclude node 10 from targets to avoid self-loop
        edges = [(10, i) for i in range(1, 10)]  # 9 edges to nodes 1-9
        edges.append((10, 11))  # Add one more to get degree 10
        taxa = {i: f"Taxon{i}" for i in range(1, 10)}  # Nodes 1-9
        taxa[11] = "Taxon10"  # Node 11 is the 10th leaf
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.degree(10) == 10
        assert 10 in net.tree_nodes
        assert 10 not in net.leaves  # Node 10 is internal
        assert len(net.leaves) == 10  # Nodes 1-9 and 11 are leaves
        # Node 10 is internal, not in taxa
        assert 10 not in net.leaves

    def test_node_degree_50(self) -> None:
        """Test node with degree 50."""
        edges = [(100, i) for i in range(1, 51)]
        taxa = {i: f"Taxon{i}" for i in range(1, 51)}
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.degree(100) == 50
        with expect_mixed_network_warning():
            assert net.validate() is True


class TestManyParallelEdges:
    """Test edge cases with many parallel edges."""

    def test_10_parallel_undirected_edges(self) -> None:
        """Test 10 parallel undirected edges between internal nodes."""
        # Parallel edges between internal nodes 1 and 2
        # Both need degree >= 3, so add additional edges
        edges = [(1, 2, i) for i in range(10)]
        edges.extend([(1, 100 + i) for i in range(3)])  # Additional edges from 1
        edges.extend([(2, 200 + i) for i in range(3)])  # Additional edges from 2
        taxa = {100 + i: f"Taxon{i}" for i in range(3)}
        taxa.update({200 + i: f"Taxon{i+3}" for i in range(3)})
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.number_of_edges() == 16  # 10 parallel + 3 + 3
        assert net.has_edge(1, 2, key=0)
        assert net.has_edge(1, 2, key=9)

    def test_10_parallel_directed_edges_to_hybrid(self) -> None:
        """Test 10 parallel directed edges to hybrid."""
        # Hybrid node 4: indegree 10, total_degree must be 11 (only 1 outgoing)
        # Node 3 needs degree >= 3
        directed_edges = [(3, 4, i) for i in range(10)]
        undirected_edges = [(4, 1), (3, 2), (3, 5)]
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa={1: "A", 2: "B", 5: "C"}
            )
        assert net.number_of_edges() == 13  # 10 directed + 3 undirected
        assert 4 in net.hybrid_nodes


class TestDeepTrees:
    """Test edge cases with deep trees."""

    def test_very_deep_tree(self) -> None:
        """Test very deep tree structure."""
        edges = []
        taxa = {}
        # Create a binary tree structure instead of a chain
        # Each internal node needs degree >= 3
        # Build a binary tree: root -> level1 -> level2 -> ... -> leaves
        root = 1
        node_id = 2
        leaves = []
        
        # Create a balanced tree with ~50 leaves
        # Each internal node needs degree >= 3, so add 3 children instead of 2
        def add_subtree(parent, depth, max_depth):
            nonlocal node_id
            if depth >= max_depth:
                # Add 3 leaves to parent to ensure degree >= 3
                for _ in range(3):
                    leaves.append(node_id)
                    edges.append((parent, node_id))
                    taxa[node_id] = f"Taxon{node_id}"
                    node_id += 1
                return
            # Add three children to ensure parent has degree >= 3
            child1 = node_id
            node_id += 1
            child2 = node_id
            node_id += 1
            child3 = node_id
            node_id += 1
            edges.append((parent, child1))
            edges.append((parent, child2))
            edges.append((parent, child3))
            add_subtree(child1, depth + 1, max_depth)
            add_subtree(child2, depth + 1, max_depth)
            add_subtree(child3, depth + 1, max_depth)
        
        # Add initial children to root (root needs degree >= 3)
        child1 = node_id
        node_id += 1
        child2 = node_id
        node_id += 1
        child3 = node_id
        node_id += 1
        edges.append((root, child1))
        edges.append((root, child2))
        edges.append((root, child3))
        add_subtree(child1, 1, 5)
        add_subtree(child2, 1, 5)
        # child3 becomes a leaf
        taxa[child3] = f"Taxon{child3}"
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        with expect_mixed_network_warning():
            assert net.validate() is True
        assert net.number_of_nodes() >= 50


class TestComplexTopologies:
    """Test edge cases with complex topologies."""

    def test_multiple_independent_hybrids(self) -> None:
        """Test multiple independent hybrid events."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Hybrid node 14: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 5, 6, 10, 15, 16, 20 need degree >= 3
        # Connect the two components
        directed_edges = [
            (5, 4), (6, 4),  # First hybrid event (corrected: 4 is hybrid)
            (15, 14), (16, 14),  # Second hybrid event (corrected: 14 is hybrid)
        ]
        undirected_edges = [
            (4, 1), (5, 8), (5, 10), (6, 9), (6, 11),  # First hybrid region
            (14, 2), (15, 18), (15, 20), (15, 22), (16, 19), (16, 21),  # Second hybrid region (15 needs degree >= 3)
            (10, 20), (10, 23), (20, 24),  # Connect the two regions, nodes 10 and 20 need degree >= 3
        ]
        taxa = {1: "A", 2: "B", 8: "C", 9: "D", 11: "E", 18: "F", 19: "G", 21: "H", 22: "I", 23: "J", 24: "K"}
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert len(net.hybrid_nodes) == 2
        with expect_mixed_network_warning():
            assert net.validate() is True

    def test_nested_hybridization(self) -> None:
        """Test nested hybridization."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Hybrid node 5: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 6, 7, 8, 9 need degree >= 3, and connect node 9
        directed_edges = [
            (7, 5), (8, 5),  # Hybrid 5
            (5, 4), (6, 4)   # Hybrid 4 (nested)
        ]
        undirected_edges = [
            (4, 1), (6, 11), (6, 12), (6, 9), (7, 13), (7, 14), (8, 15), (8, 16), (9, 2), (9, 3)
        ]
        taxa = {1: "A", 2: "B", 3: "C", 11: "D", 12: "E", 13: "F", 14: "G", 15: "H", 16: "I"}
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes
        with expect_mixed_network_warning():
            assert net.validate() is True


class TestInvalidStructures:
    """Test edge cases with invalid structures."""

    def test_disconnected_network(self) -> None:
        """Test disconnected network raises error."""
        with pytest.raises(ValueError, match="not connected"):
            MixedPhyNetwork(
                undirected_edges=[(1, 2), (3, 4)],
                taxa={2: "A", 4: "B"}
            )

    def test_internal_node_degree_two(self) -> None:
        """Test internal node with degree 2 raises error."""
        with pytest.raises(ValueError, match="degree"):
            MixedPhyNetwork(
                undirected_edges=[(1, 2), (2, 3)],
                taxa={1: "A", 3: "B"}
            )

    def test_node_invalid_indegree(self) -> None:
        """Test node with invalid indegree raises error."""
        # Node 2 has indegree 1 but total_degree 2 (should be 0 or total_degree-1=1)
        # Node 1 needs degree >= 3
        with pytest.raises(ValueError, match="indegree|degree"):
            MixedPhyNetwork(
                directed_edges=[(1, 2)],
                undirected_edges=[(2, 3), (1, 4), (1, 5)],
                taxa={3: "A", 4: "B", 5: "C"}
            )

