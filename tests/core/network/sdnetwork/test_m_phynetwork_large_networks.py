"""
Comprehensive tests for MixedPhyNetwork with large and complex networks.

This module tests large networks including:
- Binary trees with many leaves
- Networks with many hybrid nodes
- Multiple independent hybrid events
- Nested hybridization
- Stress tests
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestLargeTrees:
    """Test cases for large tree structures."""

    def test_star_tree_100_leaves(self) -> None:
        """Test star tree with 100 leaves."""
        # Create 100 leaves connected to center node 100
        # Exclude node 100 from being a target to avoid self-loop
        edges = [(100, i) for i in range(1, 100)]  # 99 edges to nodes 1-99
        # Add one more leaf to get 100 leaves total
        edges.append((100, 101))  # Node 101 is the 100th leaf
        taxa = {i: f"Taxon{i}" for i in range(1, 100)}  # Nodes 1-99
        taxa[101] = "Taxon100"  # Node 101 is the 100th leaf
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.number_of_nodes() == 101  # 100 leaves + 1 center
        assert net.number_of_edges() == 100
        assert len(net.leaves) == 100  # Nodes 1-99 and 101 are leaves
        assert 100 not in net.leaves  # Node 100 is internal (center)
        with expect_mixed_network_warning():
            assert net.validate() is True

    def test_binary_tree_structure(self) -> None:
        """Test binary tree structure."""
        edges = []
        taxa = {}
        node_counter = 1
        
        def build_tree(parent, level, max_level):
            nonlocal node_counter
            if level >= max_level:
                left_leaf = node_counter
                node_counter += 1
                right_leaf = node_counter
                node_counter += 1
                edges.append((parent, left_leaf))
                edges.append((parent, right_leaf))
                taxa[left_leaf] = f"Taxon{left_leaf}"
                taxa[right_leaf] = f"Taxon{right_leaf}"
                return
            
            left = node_counter
            node_counter += 1
            right = node_counter
            node_counter += 1
            
            edges.append((parent, left))
            edges.append((parent, right))
            build_tree(left, level + 1, max_level)
            build_tree(right, level + 1, max_level)
        
        root = 10000
        build_tree(root, 0, 5)  # 5 levels to ensure all internal nodes have degree >= 3
        # Add one more child to root to ensure it has degree >= 3
        extra_leaf = node_counter
        node_counter += 1
        edges.append((root, extra_leaf))
        taxa[extra_leaf] = f"Taxon{extra_leaf}"
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.number_of_nodes() >= 50
        with expect_mixed_network_warning():
            assert net.validate() is True


class TestManyHybrids:
    """Test cases for networks with many hybrid nodes."""

    def test_10_independent_hybrids(self) -> None:
        """Test network with 10 independent hybrid events."""
        directed_edges = []
        undirected_edges = []
        taxa = {}
        
        # Connect all components via a central node
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
        
        # Central node needs degree >= 3 (already has 20 edges from parents)
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert len(net.hybrid_nodes) == 10
        with expect_mixed_network_warning():
            assert net.validate() is True

    def test_20_hybrids(self) -> None:
        """Test network with 20 hybrid nodes."""
        directed_edges = []
        undirected_edges = []
        taxa = {}
        
        # Connect all components via a central node
        central = 0
        for i in range(20):
            parent1 = 10000 + 2 * i
            parent2 = 10000 + 2 * i + 1
            hybrid = 20000 + i
            leaf = 30000 + i
            
            directed_edges.append((parent1, hybrid))
            directed_edges.append((parent2, hybrid))
            undirected_edges.append((hybrid, leaf))
            undirected_edges.append((parent1, 40000 + 2 * i))
            undirected_edges.append((parent2, 40000 + 2 * i + 1))
            # Connect both parents to central node to ensure degree >= 3
            undirected_edges.append((central, parent1))
            undirected_edges.append((central, parent2))
            taxa[leaf] = f"Taxon{i}"
            taxa[40000 + 2 * i] = f"Taxon{20 + 2 * i}"
            taxa[40000 + 2 * i + 1] = f"Taxon{20 + 2 * i + 1}"
        
        # Central node needs degree >= 3 (already has 40 edges from parents)
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert len(net.hybrid_nodes) == 20
        with expect_mixed_network_warning():
            assert net.validate() is True


class TestNestedHybridization:
    """Test cases for nested hybridization."""

    def test_deeply_nested_hybrids(self) -> None:
        """Test deeply nested hybridization."""
        # Create chain of nested hybrids: h1 -> h2 -> h3 -> h4
        # Each hybrid: indegree 2, total_degree must be 3 (only 1 outgoing)
        # The outgoing edge from each hybrid can be either directed or undirected
        # Nodes 10, 11, 12, 13, 14 need degree >= 3
        directed_edges = [
            (10, 4), (11, 4),  # Hybrid 4: indegree 2
            (4, 3), (12, 3),   # Hybrid 3: indegree 2, 4->3 is outgoing from 4
            (3, 2), (13, 2),   # Hybrid 2: indegree 2, 3->2 is outgoing from 3
            (2, 1), (14, 1),   # Hybrid 1: indegree 2, 2->1 is outgoing from 2
        ]
        undirected_edges = [
            (1, 20),  # One outgoing from hybrid 1
            (10, 30), (10, 35), (11, 31), (11, 36),
            (12, 32), (12, 37), (13, 33), (13, 38),
            (14, 34), (14, 39)
        ]
        taxa = {20: "A",
                30: "E", 31: "F", 32: "G", 33: "H", 34: "I",
                35: "J", 36: "K", 37: "L", 38: "M", 39: "N"}
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert len(net.hybrid_nodes) == 4
        with expect_mixed_network_warning():
            assert net.validate() is True


class TestComplexTopologies:
    """Test cases for complex network topologies."""

    def test_mixed_tree_and_hybrid_structure(self) -> None:
        """Test network mixing tree and hybrid structures."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Hybrid node 5: indegree 2, total_degree must be 3 (only 1 outgoing)
        # The outgoing from 5 is (5, 4) directed, so 5 has no undirected outgoing
        # Nodes 6, 7, 10 need degree >= 3
        directed_edges = [
            (10, 5), (7, 5),  # Hybrid 5: indegree 2
            (5, 4), (6, 4),    # Hybrid 4: indegree 2, 5->4 is outgoing from 5
        ]
        undirected_edges = [
            (4, 1),  # One outgoing from hybrid 4
            (10, 11), (10, 12), (10, 13),  # Tree edges from 10
            (7, 17), (7, 18),  # Additional edges for 7
            (6, 9), (6, 16)  # Additional edges for 6
        ]
        taxa = {1: "A", 9: "B", 11: "C", 12: "D", 13: "E", 16: "F", 17: "G", 18: "H"}
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa
            )
        assert len(net.hybrid_nodes) == 2
        assert len(net.tree_nodes) >= 1
        with expect_mixed_network_warning():
            assert net.validate() is True

    def test_network_with_attributes(self) -> None:
        """Test large network with edge attributes."""
        edges = []
        taxa = {}
        for i in range(1, 51):
            edges.append({
                'u': 100,
                'v': i,
                'branch_length': i * 0.01,
                'bootstrap': 0.9 + (i % 10) * 0.01
            })
            taxa[i] = f"Taxon{i}"
        
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.number_of_nodes() == 51
        assert net.get_branch_length(100, 1) == 0.01
        with expect_mixed_network_warning():
            assert net.validate() is True

