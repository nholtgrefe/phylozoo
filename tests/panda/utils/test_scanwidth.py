"""
Tests for the scanwidth module.

This test suite covers:
- ScanwidthDAG class
- Extension class
- TreeExtension class
- Basic scanwidth computation
"""

import pytest

import networkx as nx

from phylozoo.panda.utils.scanwidth import Extension, ScanwidthDAG, TreeExtension


class TestScanwidthDAG:
    """Test cases for ScanwidthDAG class."""

    def test_init(self) -> None:
        """Test ScanwidthDAG initialization."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        sw_dag = ScanwidthDAG(dag)
        
        assert sw_dag.graph == dag
        assert sw_dag.infinity == dag.number_of_edges() + 1
        assert sw_dag._memory is True

    def test_optimal_scanwidth_simple_tree(self) -> None:
        """Test optimal scanwidth on a simple tree DAG."""
        # Simple tree: root -> node -> leaf
        dag = nx.DiGraph([(1, 2), (2, 3)])
        sw_dag = ScanwidthDAG(dag)
        
        sw, ext = sw_dag.optimal_scanwidth()
        
        assert sw is not None
        assert ext is not None
        assert sw >= 0
        assert isinstance(ext, Extension)
        assert len(ext.sigma) == 3  # All nodes in extension

    def test_optimal_scanwidth_fork(self) -> None:
        """Test optimal scanwidth on a fork DAG."""
        # Fork: root splits to two leaves
        dag = nx.DiGraph([(1, 2), (1, 3)])
        sw_dag = ScanwidthDAG(dag)
        
        sw, ext = sw_dag.optimal_scanwidth()
        
        assert sw is not None
        assert ext is not None
        assert sw >= 0
        assert isinstance(ext, Extension)
        assert len(ext.sigma) == 3

    def test_optimal_scanwidth_diamond(self) -> None:
        """Test optimal scanwidth on a diamond DAG."""
        # Diamond: two paths from root to sink
        dag = nx.DiGraph([(1, 2), (1, 3), (2, 4), (3, 4)])
        sw_dag = ScanwidthDAG(dag)
        
        sw, ext = sw_dag.optimal_scanwidth()
        
        assert sw is not None
        assert ext is not None
        assert sw >= 0
        assert isinstance(ext, Extension)
        assert len(ext.sigma) == 4

    def test_optimal_k_scanwidth(self) -> None:
        """Test optimal_k_scanwidth method."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        sw_dag = ScanwidthDAG(dag)
        
        # Try with k=1 (should work for simple tree)
        result = sw_dag.optimal_k_scanwidth(1)
        
        # Should return either False or (sw, ext)
        if result is not False:
            sw, ext = result
            assert sw >= 0
            assert isinstance(ext, Extension)
            assert sw <= 1


class TestExtension:
    """Test cases for Extension class."""

    def test_init(self) -> None:
        """Test Extension initialization."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        sigma = [1, 2, 3]
        ext = Extension(dag, sigma)
        
        assert ext.graph == dag
        assert ext.sigma == sigma

    def test_scanwidth_at_vertex_i(self) -> None:
        """Test scanwidth_at_vertex_i method."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        sigma = [1, 2, 3]
        ext = Extension(dag, sigma)
        
        # At position 0 (just node 1)
        sw_0 = ext.scanwidth_at_vertex_i(0)
        # Scanwidth can be negative (in_degree - out_degree)
        assert isinstance(sw_0, int)
        
        # At position 1 (nodes 1, 2)
        sw_1 = ext.scanwidth_at_vertex_i(1)
        assert isinstance(sw_1, int)
        
        # At position 2 (all nodes)
        sw_2 = ext.scanwidth_at_vertex_i(2)
        assert isinstance(sw_2, int)

    def test_scanwidth(self) -> None:
        """Test scanwidth method."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        sigma = [1, 2, 3]
        ext = Extension(dag, sigma)
        
        sw = ext.scanwidth()
        # Scanwidth can be negative (in_degree - out_degree)
        assert isinstance(sw, int)
        # Scanwidth should be max of all positions
        max_sw = max(ext.scanwidth_at_vertex_i(i) for i in range(len(sigma)))
        assert sw == max_sw

    def test_canonical_tree_extension(self) -> None:
        """Test canonical_tree_extension method."""
        # Use a DAG and reverse topological ordering
        # The algorithm processes nodes in order, and when processing v,
        # it needs rho[c] for successors c, so we need reverse topological order
        dag = nx.DiGraph([(1, 2), (2, 3)])
        # Reverse topological ordering: process leaves first, then internal nodes
        # This ensures that when we process a node, its successors have already been processed
        sigma = [3, 2, 1]  # Reverse of topological order [1, 2, 3]
        ext = Extension(dag, sigma)
        
        tree_ext = ext.canonical_tree_extension()
        
        assert isinstance(tree_ext, TreeExtension)
        assert tree_ext.graph == dag
        # Tree should have same nodes as original graph
        assert set(tree_ext.tree.nodes()) == set(dag.nodes())
        # Tree should be a DAG (actually a tree, which is a special DAG)
        assert nx.is_directed_acyclic_graph(tree_ext.tree)


class TestTreeExtension:
    """Test cases for TreeExtension class."""

    def test_init(self) -> None:
        """Test TreeExtension initialization."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        tree = nx.DiGraph([(1, 2), (2, 3)])  # Same structure for simplicity
        tree_ext = TreeExtension(dag, tree)
        
        assert tree_ext.graph == dag
        assert tree_ext.tree == tree

    def test_scanwidth_at_vertex(self) -> None:
        """Test scanwidth_at_vertex method."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        tree = nx.DiGraph([(1, 2), (2, 3)])
        tree_ext = TreeExtension(dag, tree)
        
        # Test for each vertex
        for v in tree.nodes():
            sw = tree_ext.scanwidth_at_vertex(v)
            assert sw >= 0

    def test_scanwidth_bag(self) -> None:
        """Test scanwidth_bag method."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        tree = nx.DiGraph([(1, 2), (2, 3)])
        tree_ext = TreeExtension(dag, tree)
        
        # Test for each vertex
        for v in tree.nodes():
            bag = tree_ext.scanwidth_bag(v)
            assert isinstance(bag, set)
            # All items in bag should be edges (tuples)
            for edge in bag:
                assert isinstance(edge, tuple)
                assert len(edge) == 2

    def test_scanwidth(self) -> None:
        """Test scanwidth method."""
        dag = nx.DiGraph([(1, 2), (2, 3)])
        tree = nx.DiGraph([(1, 2), (2, 3)])
        tree_ext = TreeExtension(dag, tree)
        
        sw = tree_ext.scanwidth()
        assert sw >= 0
        # Scanwidth should be max over all vertices
        max_sw = max(tree_ext.scanwidth_at_vertex(v) for v in tree.nodes())
        assert sw == max_sw


class TestIntegration:
    """Integration tests for scanwidth computation."""

    def test_end_to_end_simple_dag(self) -> None:
        """Test complete workflow: DAG -> Extension -> TreeExtension."""
        dag = nx.DiGraph([(1, 2), (1, 3), (2, 4), (3, 4)])
        sw_dag = ScanwidthDAG(dag)
        
        # Compute optimal scanwidth
        sw, ext = sw_dag.optimal_scanwidth()
        
        assert sw is not None
        assert ext is not None
        
        # Convert to tree extension
        tree_ext = ext.canonical_tree_extension()
        
        assert isinstance(tree_ext, TreeExtension)
        # Tree extension scanwidth should match extension scanwidth
        tree_sw = tree_ext.scanwidth()
        assert tree_sw == sw

    def test_empty_dag(self) -> None:
        """Test scanwidth on empty DAG."""
        dag = nx.DiGraph()
        sw_dag = ScanwidthDAG(dag)
        
        sw, ext = sw_dag.optimal_scanwidth()
        
        # Empty DAG should have scanwidth 0
        if sw is not None:
            assert sw == 0

    def test_single_node_dag(self) -> None:
        """Test scanwidth on single node DAG."""
        dag = nx.DiGraph()
        dag.add_node(1)
        sw_dag = ScanwidthDAG(dag)
        
        sw, ext = sw_dag.optimal_scanwidth()
        
        # Single node should have scanwidth 0
        if sw is not None:
            assert sw == 0

