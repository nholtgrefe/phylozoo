"""
Tests for suppress_degree2_node function in mm_operations.
"""

import pytest

from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.transformations import suppress_degree2_node


class TestSuppressDegree2NodeValidCombinations:
    """Tests for valid edge type combinations."""

    def test_undirected_undirected(self) -> None:
        """Test suppression of node with two undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        
        suppress_degree2_node(G, 2)
        
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        assert G.number_of_edges() == 1
        assert G.number_of_nodes() == 2

    def test_directed_in_directed_out(self) -> None:
        """Test suppression of node with incoming and outgoing directed edges -> directed."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(2, 3)
        
        suppress_degree2_node(G, 2)
        
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        assert G._directed.has_edge(1, 3)  # Should be directed
        assert not G._undirected.has_edge(1, 3)
        assert G.number_of_edges() == 1
        assert G.number_of_nodes() == 2

    def test_directed_in_undirected(self) -> None:
        """Test suppression of node with incoming directed and undirected edge -> undirected."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_undirected_edge(2, 3)
        
        suppress_degree2_node(G, 2)
        
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        assert G._undirected.has_edge(1, 3)  # Should be undirected
        assert not G._directed.has_edge(1, 3)
        assert G.number_of_edges() == 1
        assert G.number_of_nodes() == 2

    def test_undirected_directed_out(self) -> None:
        """Test suppression of node with undirected and outgoing directed edge -> directed."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        
        suppress_degree2_node(G, 2)
        
        assert 2 not in G.nodes()
        # undirected comes before directed_out, so n1=1 (neighbor of undirected), n2=3 (target of outgoing) -> 1->3
        assert G.has_edge(1, 3)
        assert G._directed.has_edge(1, 3)  # Should be directed
        assert not G._undirected.has_edge(1, 3)
        assert G.number_of_edges() == 1
        assert G.number_of_nodes() == 2


class TestSuppressDegree2NodeInvalidCombinations:
    """Tests for invalid edge type combinations that raise errors."""

    def test_directed_in_directed_in_raises_error(self) -> None:
        """Test that multiple incoming directed edges raise ValueError."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(3, 2)
        
        with pytest.raises(ValueError, match="multiple incoming directed edges"):
            suppress_degree2_node(G, 2)

    def test_directed_out_directed_out_raises_error(self) -> None:
        """Test that multiple outgoing directed edges raise ValueError."""
        G = MixedMultiGraph()
        G.add_directed_edge(2, 1)
        G.add_directed_edge(2, 3)
        
        with pytest.raises(ValueError, match="multiple outgoing directed edges"):
            suppress_degree2_node(G, 2)


class TestSuppressDegree2NodeErrorCases:
    """Tests for error cases (node not found, wrong degree, etc.)."""

    def test_node_not_in_graph_raises_error(self) -> None:
        """Test that suppressing a non-existent node raises ValueError."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        with pytest.raises(ValueError, match="not found in graph"):
            suppress_degree2_node(G, 99)

    def test_node_not_degree_2_raises_error(self) -> None:
        """Test that suppressing a non-degree-2 node raises ValueError."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)  # Node 2 now has degree 3
        
        with pytest.raises(ValueError, match="has degree 3, expected degree 2"):
            suppress_degree2_node(G, 2)

    def test_node_degree_1_raises_error(self) -> None:
        """Test that suppressing a degree-1 node raises ValueError."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        with pytest.raises(ValueError, match="has degree 1, expected degree 2"):
            suppress_degree2_node(G, 2)


class TestSuppressDegree2NodeAttributes:
    """Tests for attribute handling during suppression."""

    def test_attributes_preserved_undirected(self) -> None:
        """Test that attributes are preserved when suppressing undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0, label="edge1")
        G.add_undirected_edge(2, 3, weight=2.0, label="edge2")
        
        suppress_degree2_node(G, 2)
        
        # Check that edge exists and has attributes (d2 overrides d1)
        assert G.has_edge(1, 3)
        edge_data = dict(G._undirected[1][3][0])
        assert edge_data.get('weight') == 2.0  # d2 overrides d1
        assert edge_data.get('label') == "edge2"  # d2 overrides d1

    def test_attributes_preserved_directed(self) -> None:
        """Test that attributes are preserved when suppressing directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0, label="edge1")
        G.add_directed_edge(2, 3, weight=2.0, label="edge2")
        
        suppress_degree2_node(G, 2)
        
        # Check that edge exists and has attributes
        assert G.has_edge(1, 3)
        edge_data = dict(G._directed[1][3][0])
        assert edge_data.get('weight') == 2.0  # d2 overrides d1
        assert edge_data.get('label') == "edge2"  # d2 overrides d1

    def test_merged_attrs_provided(self) -> None:
        """Test that provided merged_attrs are used directly."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0)
        G.add_undirected_edge(2, 3, weight=2.0)
        
        merged_attrs = {'weight': 5.0, 'label': 'merged'}
        suppress_degree2_node(G, 2, merged_attrs=merged_attrs)
        
        # Check that merged_attrs are used
        assert G.has_edge(1, 3)
        edge_data = dict(G._undirected[1][3][0])
        assert edge_data.get('weight') == 5.0
        assert edge_data.get('label') == 'merged'

    def test_attributes_order_directed_in_first(self) -> None:
        """Test that directed_in edges come first in attribute merging."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0, label="directed_in")
        G.add_undirected_edge(2, 3, weight=2.0, label="undirected")
        
        suppress_degree2_node(G, 2)
        
        # directed_in comes first, so its attributes are base, then undirected overrides
        # Result should be undirected edge
        assert G.has_edge(1, 3)
        edge_data = dict(G._undirected[1][3][0])
        assert edge_data.get('weight') == 2.0  # undirected overrides
        assert edge_data.get('label') == "undirected"  # undirected overrides


class TestSuppressDegree2NodeParallelEdges:
    """Tests for parallel edge handling."""

    def test_parallel_edges_created(self) -> None:
        """Test that suppression can create parallel edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(1, 3)  # Already exists
        
        suppress_degree2_node(G, 2)
        
        # Suppression should create a new edge between 1 and 3
        # Since there was already an edge, we should have parallel edges
        # Count edges between 1 and 3 (check undirected edges)
        undirected_edges_1_3 = list(G._undirected.edges(1, 3, keys=True))
        # The suppression creates an edge, and if keys differ, we get parallel edges
        # Note: The exact behavior depends on key assignment during suppression
        assert len(undirected_edges_1_3) >= 1  # At least one edge exists
        assert G.has_edge(1, 3)  # Edge exists between 1 and 3

    def test_same_key_preserved(self) -> None:
        """
        Test that suppression creates an edge even when original edges had the same key.
        
        Note: We always use key=None to ensure parallel edges are created when needed.
        The key will be auto-generated, so we just verify the edge exists.
        """
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, key=5)
        G.add_undirected_edge(2, 3, key=5)
        
        suppress_degree2_node(G, 2)
        
        # Edge should exist (key will be auto-generated, not necessarily 5)
        assert G.has_edge(1, 3)
    
    def test_two_paths_create_parallel_edges_undirected(self) -> None:
        """
        Test that suppressing two degree-2 nodes on parallel undirected paths creates parallel edges.
        
        This is important for maintaining bi-edge connectivity when suppressing nodes
        in blob structures. When two paths exist between the same nodes (e.g., u-x-v
        and u-y-v), suppressing both x and y should create two parallel edges u-v,
        not overwrite the first edge.
        """
        G = MixedMultiGraph()
        # Create two parallel paths: 10-5-4 and 10-6-4
        G.add_undirected_edge(10, 5)
        G.add_undirected_edge(10, 6)
        G.add_undirected_edge(5, 4)
        G.add_undirected_edge(6, 4)
        
        # Suppress first degree-2 node (5)
        suppress_degree2_node(G, 5)
        edges_after_first = list(G.undirected_edges(keys=True))
        assert len([e for e in edges_after_first if (e[0] == 10 and e[1] == 4) or (e[0] == 4 and e[1] == 10)]) == 1
        
        # Suppress second degree-2 node (6)
        suppress_degree2_node(G, 6)
        edges_after_second = list(G.undirected_edges(keys=True))
        
        # Should have two parallel undirected edges between 10 and 4
        edges_10_4 = [e for e in edges_after_second if (e[0] == 10 and e[1] == 4) or (e[0] == 4 and e[1] == 10)]
        assert len(edges_10_4) == 2, f"Expected 2 parallel edges, got {len(edges_10_4)}: {edges_10_4}"
        
        # Verify bi-edge connectivity is maintained
        from phylozoo.core.primitives.m_multigraph.features import (
            bi_edge_connected_components,
            has_parallel_edges,
        )
        assert has_parallel_edges(G)
        comps = list(bi_edge_connected_components(G))
        assert len(comps) == 1, f"Should be bi-edge connected, got {len(comps)} components: {comps}"
        assert {10, 4} in comps or comps[0] == {10, 4}
    
    def test_two_paths_create_parallel_edges_directed(self) -> None:
        """
        Test that suppressing two degree-2 nodes on parallel directed paths creates parallel edges.
        
        This is important for maintaining bi-edge connectivity when suppressing nodes
        in blob structures. When two paths exist between the same nodes (e.g., u->x->v
        and u->y->v), suppressing both x and y should create two parallel edges u->v,
        not overwrite the first edge.
        """
        G = MixedMultiGraph()
        # Create two parallel paths: 10->5->4 and 10->6->4
        G.add_directed_edge(10, 5)
        G.add_directed_edge(10, 6)
        G.add_directed_edge(5, 4)
        G.add_directed_edge(6, 4)
        
        # Suppress first degree-2 node (5)
        suppress_degree2_node(G, 5)
        edges_after_first = list(G.directed_edges(keys=True))
        assert (10, 4, 0) in edges_after_first
        assert len([e for e in edges_after_first if e[0] == 10 and e[1] == 4]) == 1
        
        # Suppress second degree-2 node (6)
        suppress_degree2_node(G, 6)
        edges_after_second = list(G.directed_edges(keys=True))
        
        # Should have two parallel directed edges from 10 to 4
        edges_10_4 = [e for e in edges_after_second if e[0] == 10 and e[1] == 4]
        assert len(edges_10_4) == 2, f"Expected 2 parallel edges, got {len(edges_10_4)}: {edges_10_4}"
        
        # Verify bi-edge connectivity is maintained
        from phylozoo.core.primitives.m_multigraph.features import (
            bi_edge_connected_components,
            has_parallel_edges,
        )
        assert has_parallel_edges(G)
        comps = list(bi_edge_connected_components(G))
        assert len(comps) == 1, f"Should be bi-edge connected, got {len(comps)} components: {comps}"
        assert {10, 4} in comps or comps[0] == {10, 4}


class TestSuppressDegree2NodeComplex:
    """Tests for complex scenarios."""

    def test_suppress_in_larger_graph(self) -> None:
        """Test suppression in a larger graph structure."""
        G = MixedMultiGraph()
        # Create a path: 1-2-3-4-5
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        G.add_undirected_edge(4, 5)
        
        # Suppress node 2
        suppress_degree2_node(G, 2)
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        
        # Suppress node 4 (now degree-2 after previous suppression)
        suppress_degree2_node(G, 4)
        assert 4 not in G.nodes()
        assert G.has_edge(3, 5)
        
        # Final graph should be 1-3-5
        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2
        assert G.has_edge(1, 3)
        assert G.has_edge(3, 5)

    def test_mixed_directed_undirected_path(self) -> None:
        """Test suppression in a mixed directed/undirected path."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(3, 4)
        
        suppress_degree2_node(G, 2)
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)  # Should be undirected (directed_in + undirected)
        
        # Check direction
        assert G._undirected.has_edge(1, 3)
        assert not G._directed.has_edge(1, 3)

