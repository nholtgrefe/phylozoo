"""
Comprehensive tests for DirectedPhyNetwork with large and complex networks.

This module tests large networks including:
- Binary trees with many leaves
- Ternary trees
- Networks with many hybrid nodes
- Level-k hybridization
- Multiple independent hybrid events
- Nested hybridization
- Real-world inspired topologies
- Stress tests (1000+ nodes)
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestBinaryTrees:
    """Test cases for large binary trees."""

    def test_binary_tree_100_leaves(self) -> None:
        """Test binary tree with 100 leaves."""
        edges = []
        taxa = {}
        node_counter = 1
        
        def build_binary_tree(parent, level, max_level):
            nonlocal node_counter
            if level >= max_level:
                # At max level, create two leaves to ensure parent has 2 children
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
            build_binary_tree(left, level + 1, max_level)
            build_binary_tree(right, level + 1, max_level)
        
        root = 10000
        build_binary_tree(root, 0, 7)  # 7 levels gives up to 128 leaves
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.number_of_nodes() >= 100
        assert len(net.leaves) >= 100
        assert net.validate() is True
        assert net.is_tree() is True
        assert net.root_node == root

    def test_binary_tree_operations(self) -> None:
        """Test operations on large binary tree."""
        edges = []
        taxa = {}
        node_counter = 1
        
        def build_tree(parent, level, max_level):
            nonlocal node_counter
            if level >= max_level:
                # At max level, create two leaves to ensure parent has 2 children
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
        build_tree(root, 0, 6)  # 6 levels
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        
        # Test all operations complete
        assert net.number_of_nodes() > 0
        assert net.number_of_edges() > 0
        assert len(net.leaves) > 0
        assert len(net.taxa) > 0
        assert net.root_node == root
        assert len(net.internal_nodes) > 0
        assert net.is_tree() is True


class TestTernaryTrees:
    """Test cases for ternary trees."""

    def test_ternary_tree_100_leaves(self) -> None:
        """Test ternary tree with 100 leaves."""
        edges = []
        taxa = {}
        node_counter = 1
        
        def build_ternary_tree(parent, level, max_level):
            nonlocal node_counter
            if level >= max_level:
                # At max level, create three leaves to ensure parent has 3 children
                leaf1 = node_counter
                node_counter += 1
                leaf2 = node_counter
                node_counter += 1
                leaf3 = node_counter
                node_counter += 1
                edges.append((parent, leaf1))
                edges.append((parent, leaf2))
                edges.append((parent, leaf3))
                taxa[leaf1] = f"Taxon{leaf1}"
                taxa[leaf2] = f"Taxon{leaf2}"
                taxa[leaf3] = f"Taxon{leaf3}"
                return
            
            child1 = node_counter
            node_counter += 1
            child2 = node_counter
            node_counter += 1
            child3 = node_counter
            node_counter += 1
            
            edges.append((parent, child1))
            edges.append((parent, child2))
            edges.append((parent, child3))
            build_ternary_tree(child1, level + 1, max_level)
            build_ternary_tree(child2, level + 1, max_level)
            build_ternary_tree(child3, level + 1, max_level)
        
        root = 10000
        build_ternary_tree(root, 0, 5)  # 5 levels gives up to 243 leaves
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.number_of_nodes() >= 100
        assert len(net.leaves) >= 100
        assert net.validate() is True
        assert net.is_tree() is True


class TestNetworksWithManyHybrids:
    """Test cases for networks with many hybrid nodes."""

    def test_network_50_hybrid_nodes(self) -> None:
        """Test network with 50 hybrid nodes."""
        edges = []
        taxa = {}
        root = 10000
        
        # Create 50 independent hybrid events
        for i in range(50):
            tree1 = 1000 + 2 * i
            tree2 = 1000 + 2 * i + 1
            hybrid = 2000 + i
            leaf1 = 3000 + 2 * i
            leaf2 = 3000 + 2 * i + 1
            
            edges.append((root, tree1))
            edges.append((root, tree2))
            edges.append((tree1, hybrid))
            edges.append((tree1, leaf1))  # tree1 splits to hybrid and leaf1
            edges.append((tree2, hybrid))
            edges.append((tree2, leaf2))  # tree2 splits to hybrid and leaf2
            edges.append((hybrid, 4000 + i))  # hybrid has one child
            taxa[leaf1] = f"Taxon{2*i}"
            taxa[leaf2] = f"Taxon{2*i+1}"
            taxa[4000 + i] = f"Taxon{100+i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert len(net.hybrid_nodes) == 50
        assert net.validate() is True
        assert net.is_tree() is False

    def test_network_100_hybrid_nodes(self) -> None:
        """Test network with 100 hybrid nodes."""
        edges = []
        taxa = {}
        root = 100000
        
        for i in range(100):
            tree1 = 10000 + 2 * i
            tree2 = 10000 + 2 * i + 1
            hybrid = 20000 + i
            leaf1 = 30000 + 2 * i
            leaf2 = 30000 + 2 * i + 1
            
            edges.append((root, tree1))
            edges.append((root, tree2))
            edges.append((tree1, hybrid))
            edges.append((tree1, leaf1))  # tree1 splits to hybrid and leaf1
            edges.append((tree2, hybrid))
            edges.append((tree2, leaf2))  # tree2 splits to hybrid and leaf2
            edges.append((hybrid, 40000 + i))  # hybrid has one child
            taxa[leaf1] = f"Taxon{2*i}"
            taxa[leaf2] = f"Taxon{2*i+1}"
            taxa[40000 + i] = f"TaxonH{i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert len(net.hybrid_nodes) == 100
        assert net.validate() is True


class TestLevelKHybridization:
    """Test cases for level-k hybridization."""

    def test_level_2_hybridization(self) -> None:
        """Test network with level-2 hybridization."""
        # Two hybrid nodes that share a parent
        edges = [
            (10, 7), (10, 8),  # Root splits
            (7, 5), (7, 6),    # Tree nodes
            (8, 5), (8, 9),    # Hybrid 5
            (5, 4), (6, 4),    # Hybrid 4 (child of hybrid 5)
            (6, 11),  # 6 splits to 4 and 11
            (4, 1), (9, 2), (9, 3)  # To leaves (9 splits to 2 children)
        ]
        net = DirectedPhyNetwork(edges=edges, taxa={1: "A", 2: "B", 3: "C", 11: "D"})
        assert len(net.hybrid_nodes) == 2
        assert net.validate() is True

    def test_level_3_hybridization(self) -> None:
        """Test network with level-3 hybridization."""
        # Three levels of hybridization
        edges = [
            (20, 15), (20, 16),  # Root
            (15, 10), (15, 11),  # Tree nodes
            (16, 10), (16, 12),  # Hybrid 10
            (10, 5), (11, 5),    # Hybrid 5 (child of hybrid 10)
            (5, 4), (12, 4),     # Hybrid 4 (child of hybrid 5)
            (11, 13),  # 11 splits to 5 and 13
            (12, 14),  # 12 splits to 4 and 14
            (13, 2), (13, 17),  # 13 splits to 2 children
            (14, 3), (14, 18),  # 14 splits to 2 children
            (4, 1)  # 4 is hybrid, so it has one child
        ]
        net = DirectedPhyNetwork(edges=edges, taxa={1: "A", 2: "B", 3: "C", 17: "D", 18: "E"})
        assert len(net.hybrid_nodes) == 3
        assert net.validate() is True

    def test_level_4_hybridization(self) -> None:
        """Test network with level-4 hybridization."""
        # Four levels of hybridization
        edges = [
            (30, 25), (30, 26),
            (25, 20), (25, 21),
            (26, 20), (26, 22),  # Hybrid 20
            (20, 15), (21, 15),  # Hybrid 15
            (21, 23),  # 21 splits to 15 and 23
            (15, 10), (22, 10),  # Hybrid 10
            (22, 24),  # 22 splits to 10, 5, 24
            (10, 5), (22, 5),    # Hybrid 5
            (5, 1), (23, 2), (23, 4), (24, 3), (24, 6)  # To leaves (23 and 24 split to 2 children)
        ]
        net = DirectedPhyNetwork(edges=edges, taxa={1: "A", 2: "B", 3: "C", 4: "D", 6: "E"})
        assert len(net.hybrid_nodes) == 4
        assert net.validate() is True


class TestMultipleIndependentHybrids:
    """Test cases for multiple independent hybrid events."""

    def test_10_independent_hybrids(self) -> None:
        """Test 10 independent hybrid events."""
        edges = []
        taxa = {}
        root = 1000
        
        for i in range(10):
            t1 = 100 + 2 * i
            t2 = 100 + 2 * i + 1
            h = 200 + i
            l1 = 300 + 2 * i
            l2 = 300 + 2 * i + 1
            
            edges.append((root, t1))
            edges.append((root, t2))
            edges.append((t1, h))
            edges.append((t1, l1))  # t1 splits to h and l1
            edges.append((t2, h))
            edges.append((t2, l2))  # t2 splits to h and l2
            edges.append((h, 400 + i))  # h has one child
            taxa[l1] = f"Taxon{2*i}"
            taxa[l2] = f"Taxon{2*i+1}"
            taxa[400 + i] = f"TaxonH{i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert len(net.hybrid_nodes) == 10
        assert net.validate() is True

    def test_20_independent_hybrids(self) -> None:
        """Test 20 independent hybrid events."""
        edges = []
        taxa = {}
        root = 10000
        
        for i in range(20):
            t1 = 1000 + 2 * i
            t2 = 1000 + 2 * i + 1
            h = 2000 + i
            l1 = 3000 + 2 * i
            l2 = 3000 + 2 * i + 1
            
            edges.append((root, t1))
            edges.append((root, t2))
            edges.append((t1, h))
            edges.append((t1, l1))  # t1 splits to h and l1
            edges.append((t2, h))
            edges.append((t2, l2))  # t2 splits to h and l2
            edges.append((h, 4000 + i))  # h has one child
            taxa[l1] = f"Taxon{2*i}"
            taxa[l2] = f"Taxon{2*i+1}"
            taxa[4000 + i] = f"TaxonH{i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert len(net.hybrid_nodes) == 20
        assert net.validate() is True


class TestNestedHybridization:
    """Test cases for nested hybridization."""

    def test_deeply_nested_hybrids(self) -> None:
        """Test deeply nested hybridization."""
        # Hybrid 4 is child of hybrid 5, which is child of hybrid 6
        edges = [
            (20, 15), (20, 16),
            (15, 10), (15, 11),
            (16, 10), (16, 12),  # Hybrid 10
            (10, 5), (11, 5),    # Hybrid 5
            (11, 13),  # 11 splits to 5 and 13
            (5, 4), (12, 4),     # Hybrid 4
            (12, 14),  # 12 splits to 4 and 14
            (4, 1), (13, 2), (13, 3), (14, 6), (14, 7)  # To leaves (4 is hybrid, so it has one child)
        ]
        net = DirectedPhyNetwork(edges=edges, taxa={1: "A", 2: "B", 3: "C", 6: "D", 7: "E"})
        assert len(net.hybrid_nodes) == 3
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes
        assert 10 in net.hybrid_nodes
        assert net.validate() is True


class TestRealWorldInspiredTopologies:
    """Test cases inspired by real-world phylogenetic networks."""

    def test_caterpillar_tree(self) -> None:
        """Test caterpillar tree topology (one long branch with leaves)."""
        # Structure: root -> chain -> leaves (each internal node splits to next node and a leaf)
        edges = []
        taxa = {}
        root = 100
        current = root
        
        for i in range(1, 51):
            next_node = 100 + i
            edges.append((current, next_node))
            edges.append((current, 200 + i))  # Leaf
            taxa[200 + i] = f"Taxon{i}"
            current = next_node
        
        # Last node splits to 2 leaves
        edges.append((current, 200 + 51))
        edges.append((current, 200 + 52))
        taxa[200 + 51] = "Taxon51"
        taxa[200 + 52] = "Taxon52"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.validate() is True
        assert len(net.leaves) == 52

    def test_balanced_tree(self) -> None:
        """Test perfectly balanced binary tree."""
        edges = []
        taxa = {}
        node_counter = 1
        
        def build_balanced(parent, level, max_level):
            nonlocal node_counter
            if level >= max_level:
                # At max level, create two leaves to ensure parent has 2 children
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
            build_balanced(left, level + 1, max_level)
            build_balanced(right, level + 1, max_level)
        
        root = 10000
        build_balanced(root, 0, 7)  # 7 levels = 128 leaves
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.validate() is True
        assert net.is_tree() is True

    def test_network_with_mixed_topology(self) -> None:
        """Test network with mixed tree and hybrid regions."""
        edges = []
        taxa = {}
        root = 1000
        
        # Tree region
        edges.append((root, 100))
        edges.append((root, 101))
        edges.append((100, 1))
        edges.append((100, 2))
        edges.append((101, 3))
        edges.append((101, 4))
        taxa.update({1: "A", 2: "B", 3: "C", 4: "D"})
        
        # Hybrid region
        edges.append((root, 200))
        edges.append((root, 201))
        edges.append((200, 300))
        edges.append((200, 6))  # 200 splits to 300 and 6
        edges.append((201, 300))
        edges.append((201, 7))  # 201 splits to 300 and 7
        edges.append((300, 5))
        taxa.update({5: "E", 6: "F", 7: "G"})
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.validate() is True
        assert len(net.hybrid_nodes) == 1
        assert len(net.tree_nodes) > 0


class TestStressTests:
    """Stress tests with very large networks."""

    def test_1000_nodes(self) -> None:
        """Test network with 1000 nodes."""
        # Create a wide tree with 1000 nodes
        edges = []
        taxa = {}
        root = 10000
        
        # Root with many children, each leading to a small subtree
        for i in range(100):
            child = 20000 + i
            edges.append((root, child))
            # Each child has 9 leaves
            for j in range(9):
                leaf = 30000 + i * 9 + j
                edges.append((child, leaf))
                taxa[leaf] = f"Taxon{i}_{j}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.number_of_nodes() == 1001  # 1 root + 100 internal + 900 leaves
        assert net.validate() is True

    def test_2000_edges(self) -> None:
        """Test network with 2000 edges."""
        edges = []
        taxa = {}
        root = 100000
        
        # Create structure with many edges
        node_id = 1
        for level in range(10):
            for i in range(2 ** level):
                if level < 9:
                    parent = 100000 + level * 1000 + i
                    child1 = 100000 + (level + 1) * 1000 + 2 * i
                    child2 = 100000 + (level + 1) * 1000 + 2 * i + 1
                    edges.append((parent, child1))
                    edges.append((parent, child2))
                else:
                    parent = 100000 + level * 1000 + i
                    leaf1 = node_id
                    node_id += 1
                    leaf2 = node_id
                    node_id += 1
                    edges.append((parent, leaf1))
                    edges.append((parent, leaf2))
                    taxa[leaf1] = f"Taxon{leaf1}"
                    taxa[leaf2] = f"Taxon{leaf2}"
                    if len(edges) >= 2000:
                        break
            if len(edges) >= 2000:
                break
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.number_of_edges() >= 2000
        assert net.validate() is True

    def test_large_network_all_operations(self) -> None:
        """Test all operations on large network."""
        # Create large network
        edges = []
        taxa = {}
        root = 100000
        
        # Create 500 independent hybrid events
        for i in range(500):
            t1 = 200000 + 2 * i
            t2 = 200000 + 2 * i + 1
            h = 300000 + i
            l1 = 400000 + 2 * i
            l2 = 400000 + 2 * i + 1
            
            edges.append((root, t1))
            edges.append((root, t2))
            edges.append((t1, h))
            edges.append((t1, l1))  # t1 splits to h and l1
            edges.append((t2, h))
            edges.append((t2, l2))  # t2 splits to h and l2
            edges.append((h, 500000 + i))  # h has one child
            taxa[l1] = f"Taxon{2*i}"
            taxa[l2] = f"Taxon{2*i+1}"
            taxa[500000 + i] = f"TaxonH{i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        
        # Test all operations complete (may be slow, but should complete)
        assert net.number_of_nodes() > 0
        assert net.number_of_edges() > 0
        assert len(net.leaves) == 1500  # 1000 from t1/t2 + 500 from h
        assert len(net.taxa) == 1500
        assert net.root_node == root
        assert len(net.hybrid_nodes) == 500
        assert len(net.tree_nodes) == 1000
        assert len(net.hybrid_edges) == 1000
        assert net.is_tree() is False
        assert net.validate() is True
        
        # Test iteration
        nodes = list(net)
        assert len(nodes) == net.number_of_nodes()
        
        # Test degree operations
        for node in list(net._graph.nodes)[:10]:  # Sample first 10
            _ = net.degree(node)
            _ = net.indegree(node)
            _ = net.outdegree(node)


class TestComplexRealWorldScenarios:
    """Test complex scenarios inspired by real phylogenetic networks."""

    def test_network_with_multiple_hybrid_levels(self) -> None:
        """Test network with hybrid nodes at multiple levels."""
        edges = []
        taxa = {}
        root = 10000
        
        # Level 1: 5 hybrids
        for i in range(5):
            t1 = 1000 + 2 * i
            t2 = 1000 + 2 * i + 1
            h1 = 2000 + i
            l1 = 3000 + 2 * i
            l1b = 3000 + 2 * i + 1
            
            edges.append((root, t1))
            edges.append((root, t2))
            edges.append((t1, h1))
            edges.append((t1, l1))  # t1 splits to h1 and l1
            edges.append((t2, h1))
            edges.append((t2, l1b))  # t2 splits to h1 and l1b
            edges.append((h1, 7000 + i))  # h1 has one child (not a leaf, will have children)
            taxa[l1] = f"Taxon{2*i}"
            taxa[l1b] = f"Taxon{2*i+1}"
        
        # Level 2: 3 hybrids (children of level 1 hybrid children)
        for i in range(3):
            h1_child = 7000 + i  # Child of h1
            t3 = 4000 + 2 * i
            t4 = 4000 + 2 * i + 1
            h2 = 5000 + i
            l2 = 6000 + 2 * i
            l2b = 6000 + 2 * i + 1
            
            edges.append((h1_child, t3))
            edges.append((h1_child, t4))
            edges.append((t3, h2))
            edges.append((t3, l2))  # t3 splits to h2 and l2
            edges.append((t4, h2))
            edges.append((t4, l2b))  # t4 splits to h2 and l2b
            edges.append((h2, 8000 + i))  # h2 has one child (leaf)
            taxa[l2] = f"Taxon2_{2*i}"
            taxa[l2b] = f"Taxon2_{2*i+1}"
            taxa[8000 + i] = f"TaxonH2{i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert len(net.hybrid_nodes) == 8  # 5 level-1 + 3 level-2
        assert net.validate() is True

    def test_network_with_hybrid_and_tree_regions(self) -> None:
        """Test network with both hybrid and tree regions."""
        edges = []
        taxa = {}
        root = 1000
        
        # Tree region: balanced binary tree
        def add_tree(parent, level, max_level, node_id):
            if level >= max_level:
                # At max level, create two leaves to ensure parent has out-degree >= 2
                leaf1 = node_id[0]
                node_id[0] += 1
                leaf2 = node_id[0]
                node_id[0] += 1
                edges.append((parent, leaf1))
                edges.append((parent, leaf2))
                taxa[leaf1] = f"TreeTaxon{leaf1}"
                taxa[leaf2] = f"TreeTaxon{leaf2}"
                return
            
            left = node_id[0]
            node_id[0] += 1
            right = node_id[0]
            node_id[0] += 1
            
            edges.append((parent, left))
            edges.append((parent, right))
            add_tree(left, level + 1, max_level, node_id)
            add_tree(right, level + 1, max_level, node_id)
        
        tree_node_id = [1]
        add_tree(root, 0, 4, tree_node_id)  # Tree with leaves 1-15
        
        # Hybrid region: multiple hybrids
        for i in range(5):
            t1 = 100 + 2 * i
            t2 = 100 + 2 * i + 1
            h = 200 + i
            l1 = 300 + 2 * i
            l2 = 300 + 2 * i + 1
            
            edges.append((root, t1))
            edges.append((root, t2))
            edges.append((t1, h))
            edges.append((t1, l1))  # t1 splits to h and l1
            edges.append((t2, h))
            edges.append((t2, l2))  # t2 splits to h and l2
            edges.append((h, 400 + i))  # h has one child
            taxa[l1] = f"HybridTaxon{2*i}"
            taxa[l2] = f"HybridTaxon{2*i+1}"
            taxa[400 + i] = f"HybridTaxonH{i}"
        
        net = DirectedPhyNetwork(edges=edges, taxa=taxa)
        assert net.validate() is True
        assert len(net.hybrid_nodes) == 5
        assert len(net.tree_nodes) > 0

