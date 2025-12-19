"""
Tests for bidirectional directed edge prevention in identify_vertices for MixedMultiGraph.
"""

import pytest

from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.transformations import identify_vertices


class TestIdentifyVerticesBidirectionalDirectedError:
    """Tests for bidirectional directed edge error cases in MixedMultiGraph."""

    def test_bidirectional_directed_error_simple(self) -> None:
        """Test that creating bidirectional directed edges raises error."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 3)  # 1 -> 3
        G.add_directed_edge(2, 3)  # 2 -> 3
        G.add_directed_edge(3, 1)  # 3 -> 1
        
        # After merging 1 and 2: would have 1 -> 3 and 3 -> 1 (bidirectional!)
        with pytest.raises(ValueError, match="edges in both directions"):
            identify_vertices(G, [1, 2])

    def test_bidirectional_directed_error_with_existing_edge(self) -> None:
        """Test bidirectional error when first_vertex already has an edge."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 4)  # 1 -> 4
        G.add_directed_edge(2, 5)  # 2 -> 5
        G.add_directed_edge(5, 1)  # 5 -> 1
        
        # After merging 1 and 2: would have 1 -> 4, 1 -> 5, and 5 -> 1 (bidirectional!)
        with pytest.raises(ValueError, match="edges in both directions"):
            identify_vertices(G, [1, 2])

    def test_bidirectional_directed_error_complex(self) -> None:
        """Test bidirectional error with multiple vertices."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 4)  # 1 -> 4
        G.add_directed_edge(2, 5)  # 2 -> 5
        G.add_directed_edge(3, 4)  # 3 -> 4
        G.add_directed_edge(5, 1)  # 5 -> 1
        
        # After merging [1, 2, 3]: would have 1 -> 4, 1 -> 5, and 5 -> 1 (bidirectional!)
        with pytest.raises(ValueError, match="edges in both directions"):
            identify_vertices(G, [1, 2, 3])

    def test_bidirectional_directed_error_reverse(self) -> None:
        """Test bidirectional error in reverse direction."""
        G = MixedMultiGraph()
        G.add_directed_edge(4, 1)  # 4 -> 1
        G.add_directed_edge(2, 4)  # 2 -> 4
        
        # After merging 1 and 2: would have 4 -> 1 and 1 -> 4 (bidirectional!)
        with pytest.raises(ValueError, match="edges in both directions"):
            identify_vertices(G, [1, 2])

    def test_no_error_when_no_bidirectional(self) -> None:
        """Test that unidirectional edges don't raise error."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 3)  # 1 -> 3
        G.add_directed_edge(2, 3)  # 2 -> 3
        G.add_directed_edge(2, 4)  # 2 -> 4
        
        # After merging 1 and 2: would have 1 -> 3, 1 -> 4 (no bidirectional)
        identify_vertices(G, [1, 2])
        
        assert 1 in G.nodes()
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        assert G.has_edge(1, 4)

