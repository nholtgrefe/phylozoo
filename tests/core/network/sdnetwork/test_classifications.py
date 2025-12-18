"""
Tests for classification functions for semi-directed and mixed phylogenetic networks.

This module tests classification functions including:
- level
- vertex_level
- reticulation_number
- is_binary
- is_tree
"""

import warnings

import pytest

from phylozoo.core.network import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import (
    is_binary,
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
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert level(net) == 0

    def test_level_tree(self) -> None:
        """Test level in tree (no hybrids)."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert level(net) == 0

    def test_level_single_hybrid(self) -> None:
        """Test level with single hybrid node."""
        # Network with hybrid node 5
        # Hybrid edges: (6, 5), (7, 5)
        # Blob contains {5, 6, 7}
        # Hybrid edges in blob: 2
        # Hybrid nodes in blob: 1
        # Level = 2 - 1 = 1
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ]
        )
        assert level(net) == 1

    def test_level_multiple_blobs(self) -> None:
        """Test level with multiple blobs."""
        # Two separate hybrid nodes in different blobs, connected via bridge edges
        # With bi-edge connectivity, bridge edges separate the blobs
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5},  # Hybrid node 4
                {'u': 8, 'v': 7, 'gamma': 0.5}, {'u': 9, 'v': 7, 'gamma': 0.5},  # Hybrid node 7
            ],
            undirected_edges=[
                (4, 1), (5, 2), (5, 6), (6, 3),  # Blob around hybrid 4
                (7, 10), (8, 11), (8, 12), (9, 13), (9, 12),  # Blob around hybrid 7
                (6, 14), (12, 14), (14, 15),  # Bridge edges connecting via node 14
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'}), (13, {'label': 'F'}),
                (14, {'label': 'G'}), (15, {'label': 'H'}),
            ]
        )
        # With bi-edge connectivity, edges (6,14) and (12,14) are bridges,
        # so we get separate blobs. Each blob should have level 1 (2 edges - 1 node)
        assert level(net) == 1


class TestVertexLevel:
    """Test cases for vertex_level() function."""

    def test_vertex_level_empty_network(self) -> None:
        """Test vertex_level in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert vertex_level(net) == 0

    def test_vertex_level_tree(self) -> None:
        """Test vertex_level in tree (no hybrids)."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert vertex_level(net) == 0

    def test_vertex_level_single_hybrid(self) -> None:
        """Test vertex_level with single hybrid node."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ]
        )
        # Blob contains hybrid node 5
        assert vertex_level(net) == 1

    def test_vertex_level_multiple_hybrids_same_blob(self) -> None:
        """Test vertex_level with multiple hybrids in same blob."""
        # Network with nested hybrids in same blob - need to connect them properly
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 7, 'v': 5, 'gamma': 0.5}, {'u': 8, 'v': 5, 'gamma': 0.5},  # Hybrid node 5
                {'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5},  # Hybrid node 4
            ],
            undirected_edges=[
                (4, 1), (6, 2), (6, 3), (6, 9),  # Connect 6 to blob
                (7, 10), (7, 11), (8, 12), (8, 13),  # Connect 7, 8
                (9, 10), (9, 14), (10, 15), (11, 12), (11, 16), (12, 18),  # Connect components to form single blob
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (9, {'label': 'D'}), (10, {'label': 'E'}), (11, {'label': 'F'}),
                (12, {'label': 'G'}), (13, {'label': 'H'}), (14, {'label': 'I'}),
                (15, {'label': 'J'}), (16, {'label': 'K'}), (18, {'label': 'M'}),
            ]
        )
        # Blob contains hybrid nodes 4 and 5
        assert vertex_level(net) == 2


class TestReticulationNumber:
    """Test cases for reticulation_number() function."""

    def test_reticulation_number_empty_network(self) -> None:
        """Test reticulation_number in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert reticulation_number(net) == 0

    def test_reticulation_number_tree(self) -> None:
        """Test reticulation_number in tree (no hybrids)."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert reticulation_number(net) == 0

    def test_reticulation_number_single_hybrid(self) -> None:
        """Test reticulation_number with single hybrid node."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5}
            ],
            undirected_edges=[(4, 1), (5, 2), (5, 7), (6, 3), (6, 7), (7, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                   (7, {'label': 'D'}), (8, {'label': 'E'})]
        )
        # 2 hybrid edges - 1 hybrid node = 1
        assert reticulation_number(net) == 1

    def test_reticulation_number_multiple_hybrids(self) -> None:
        """Test reticulation_number with multiple hybrid nodes."""
        # Two separate hybrid nodes, not nested
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5},  # Hybrid node 4
                {'u': 8, 'v': 7, 'gamma': 0.5}, {'u': 9, 'v': 7, 'gamma': 0.5},  # Hybrid node 7
            ],
            undirected_edges=[
                (4, 1), (5, 2), (5, 6), (6, 3),  # Blob around hybrid 4
                (7, 10), (8, 11), (8, 12), (9, 13), (9, 12),  # Blob around hybrid 7
                (6, 14), (12, 14), (14, 15),  # Connect via bridge node 14
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'}), (13, {'label': 'F'}),
                (14, {'label': 'G'}), (15, {'label': 'H'}),
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
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert is_binary(net) is True

    def test_is_binary_simple_tree(self) -> None:
        """Test is_binary in simple binary tree."""
        # Internal node 3 has degree 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert is_binary(net) is True

    def test_is_binary_non_binary_internal(self) -> None:
        """Test is_binary with internal node having wrong degree."""
        # Internal node 3 has degree 4, should be 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4), (3, 5)],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (4, {'label': 'C'}), (5, {'label': 'D'})
            ]
        )
        assert is_binary(net) is False

    def test_is_binary_with_hybrids(self) -> None:
        """Test is_binary in network with hybrid nodes."""
        # Hybrid node 4 has degree 3 (in-degree 2, undirected degree 1)
        # Tree nodes 5, 6 have degree 3
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5}
            ],
            undirected_edges=[(4, 1), (5, 2), (5, 7), (6, 3), (6, 7), (7, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                   (7, {'label': 'D'}), (8, {'label': 'E'})]
        )
        assert is_binary(net) is True

    def test_is_binary_single_node(self) -> None:
        """Test is_binary in single-node network."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[], undirected_edges=[],
            nodes=[(1, {'label': 'A'})]
        )
        assert is_binary(net) is True


class TestIsTree:
    """Test cases for is_tree() function."""

    def test_is_tree_empty_network(self) -> None:
        """Test is_tree in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert is_tree(net) is True

    def test_is_tree_simple_tree(self) -> None:
        """Test is_tree in simple tree."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert is_tree(net) is True

    def test_is_tree_with_hybrids(self) -> None:
        """Test is_tree in network with hybrid nodes."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5}
            ],
            undirected_edges=[(4, 1), (5, 2), (5, 7), (6, 3), (6, 7), (7, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                   (7, {'label': 'D'}), (8, {'label': 'E'})]
        )
        assert is_tree(net) is False

    def test_is_tree_large_tree(self) -> None:
        """Test is_tree in large binary tree."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                (7, 5), (7, 6), (7, 8),
                (5, 3), (5, 4),
                (6, 1), (6, 2)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'}), (8, {'label': 'E'})
            ]
        )
        assert is_tree(net) is True

