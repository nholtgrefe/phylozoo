"""
Reusable network structures for testing and examples.

This module provides a collection of pre-defined network structures that can be
reused across tests, examples, and demonstrations. This avoids recreating the
same network structures repeatedly.

All networks are valid according to DirectedPhyNetwork, MixedPhyNetwork, and
SemiDirectedPhyNetwork validation rules.
"""

from __future__ import annotations

from typing import Any

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork


# ========== Minimal Networks ==========


def empty_network() -> DirectedPhyNetwork:
    """Empty network (no nodes, no edges)."""
    return DirectedPhyNetwork(edges=[], nodes=[])


def single_node_network() -> DirectedPhyNetwork:
    """Single-node network (root and leaf are the same)."""
    return DirectedPhyNetwork(nodes=[(1, {'label': 'A'})])


def single_edge_network() -> DirectedPhyNetwork:
    """Network with single edge (root -> leaf)."""
    return DirectedPhyNetwork(
        edges=[(1, 2)],
        nodes=[(2, {'label': 'A'})]
    )


# ========== Simple Trees ==========


def two_leaf_star() -> DirectedPhyNetwork:
    """Star tree with root and two leaves."""
    return DirectedPhyNetwork(
        edges=[(3, 1), (3, 2)],
        nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    )


def three_leaf_star() -> DirectedPhyNetwork:
    """Star tree with root and three leaves."""
    return DirectedPhyNetwork(
        edges=[(4, 1), (4, 2), (4, 3)],
        nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
    )


def binary_tree_4_leaves() -> DirectedPhyNetwork:
    """Balanced binary tree with 4 leaves."""
    return DirectedPhyNetwork(
        edges=[
            (7, 5), (7, 6),
            (5, 1), (5, 2),
            (6, 3), (6, 4),
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (4, {'label': 'D'}),
        ]
    )


def binary_tree_with_internal() -> DirectedPhyNetwork:
    """Binary tree with internal node."""
    return DirectedPhyNetwork(
        edges=[
            (4, 3), (3, 1), (3, 2), (4, 5)
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (5, {'label': 'C'}),
        ]
    )


# ========== Networks with Hybrid Nodes ==========


def single_hybrid_network() -> DirectedPhyNetwork:
    """Network with single hybrid node."""
    return DirectedPhyNetwork(
        edges=[
            (8, 5), (8, 6),
            (5, 1), (5, 2),
            (6, 3), (6, 9),
            {'u': 5, 'v': 4, 'gamma': 0.6},
            {'u': 6, 'v': 4, 'gamma': 0.4},
            (4, 10),
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (9, {'label': 'D'}),
            (10, {'label': 'E'}),
        ]
    )


def single_hybrid_network_simple() -> DirectedPhyNetwork:
    """Simpler network with single hybrid node."""
    return DirectedPhyNetwork(
        edges=[
            (7, 5), (7, 6),
            (5, 4), (5, 8),
            (6, 4), (6, 9),
            (4, 2),
        ],
        nodes=[
            (2, {'label': 'A'}),
            (8, {'label': 'B'}),
            (9, {'label': 'C'}),
        ]
    )


def two_hybrid_network() -> DirectedPhyNetwork:
    """Network with two hybrid nodes."""
    return DirectedPhyNetwork(
        edges=[
            (10, 8), (10, 9),
            (8, 5), (8, 6),
            (9, 7), (9, 14),
            (5, 1), (5, 2),
            (6, 3),
            (7, 15), (7, 16),
            {'u': 5, 'v': 4, 'gamma': 0.6},
            {'u': 6, 'v': 4, 'gamma': 0.4},
            (4, 11),
            (11, 12), (11, 18), (11, 19),
            (12, 21),
            (18, 20),
            {'u': 12, 'v': 13, 'gamma': 0.7},
            {'u': 18, 'v': 13, 'gamma': 0.3},
            (13, 17),
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (15, {'label': 'D'}),
            (16, {'label': 'E'}),
            (17, {'label': 'F'}),
            (19, {'label': 'G'}),
            (20, {'label': 'H'}),
            (21, {'label': 'I'}),
        ]
    )


def hybrid_with_parallel_edges() -> DirectedPhyNetwork:
    """Network with hybrid node and parallel edges."""
    return DirectedPhyNetwork(
        edges=[
            (7, 5), (7, 6),
            (5, 4, 0), (5, 4, 1),  # Parallel edges to hybrid
            (5, 8),
            (6, 4),
            (6, 9),
            (4, 2),
        ],
        nodes=[
            (2, {'label': 'A'}),
            (8, {'label': 'B'}),
            (9, {'label': 'C'}),
        ]
    )


# ========== Networks with Branch Lengths ==========


def tree_with_branch_lengths() -> DirectedPhyNetwork:
    """Tree with branch lengths on all edges."""
    return DirectedPhyNetwork(
        edges=[
            {'u': 7, 'v': 5, 'branch_length': 0.5},
            {'u': 7, 'v': 6, 'branch_length': 0.3},
            {'u': 5, 'v': 1, 'branch_length': 0.2},
            {'u': 5, 'v': 2, 'branch_length': 0.4},
            {'u': 6, 'v': 3, 'branch_length': 0.3},
            {'u': 6, 'v': 4, 'branch_length': 0.25},
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (4, {'label': 'D'}),
        ]
    )


def network_with_bootstrap() -> DirectedPhyNetwork:
    """Network with bootstrap support values."""
    return DirectedPhyNetwork(
        edges=[
            {'u': 7, 'v': 5, 'bootstrap': 0.95},
            {'u': 7, 'v': 6, 'bootstrap': 0.87},
            {'u': 5, 'v': 1, 'bootstrap': 0.92},
            {'u': 5, 'v': 2, 'bootstrap': 0.88},
            {'u': 6, 'v': 3, 'bootstrap': 0.91},
            {'u': 6, 'v': 4, 'bootstrap': 0.85},
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (4, {'label': 'D'}),
        ]
    )


# ========== Semi-Directed Networks ==========


def semidirected_simple() -> SemiDirectedPhyNetwork:
    """Simple semi-directed network (tree-like)."""
    return SemiDirectedPhyNetwork(
        undirected_edges=[
            (3, 1), (3, 2), (3, 4),
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (4, {'label': 'C'}),
        ]
    )


def semidirected_with_hybrid() -> SemiDirectedPhyNetwork:
    """Semi-directed network with hybrid node."""
    return SemiDirectedPhyNetwork(
        undirected_edges=[
            (5, 1),
            (6, 2), (6, 3), (6, 4),
            (6, 7), (7, 8),
        ],
        directed_edges=[
            {'u': 6, 'v': 5, 'gamma': 0.6},
            {'u': 7, 'v': 5, 'gamma': 0.4},
        ],
        nodes=[
            (1, {'label': 'A'}),
            (2, {'label': 'B'}),
            (3, {'label': 'C'}),
            (4, {'label': 'D'}),
            (8, {'label': 'E'}),
        ]
    )


# ========== Large Networks ==========


def star_tree_10_leaves() -> DirectedPhyNetwork:
    """Star tree with 10 leaves."""
    edges = [(100, i) for i in range(1, 11)]
    nodes = [(i, {'label': f'Taxon{i}'}) for i in range(1, 11)]
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


def star_tree_50_leaves() -> DirectedPhyNetwork:
    """Star tree with 50 leaves."""
    edges = [(1000, i) for i in range(1, 51)]
    nodes = [(i, {'label': f'Taxon{i}'}) for i in range(1, 51)]
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


def balanced_binary_tree(depth: int = 6) -> DirectedPhyNetwork:
    """
    Generate a balanced binary tree with specified depth.
    
    Parameters
    ----------
    depth : int, optional
        Depth of the tree. By default 6 (gives ~64 leaves).
    
    Returns
    -------
    DirectedPhyNetwork
        Balanced binary tree.
    """
    edges = []
    nodes = []
    node_counter = 1
    
    def build_tree(parent: int, level: int, max_level: int) -> None:
        nonlocal node_counter
        if level >= max_level:
            # At max level, create two leaves
            left_leaf = node_counter
            node_counter += 1
            right_leaf = node_counter
            node_counter += 1
            edges.append((parent, left_leaf))
            edges.append((parent, right_leaf))
            nodes.append((left_leaf, {'label': f'Taxon{left_leaf}'}))
            nodes.append((right_leaf, {'label': f'Taxon{right_leaf}'}))
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
    build_tree(root, 0, depth)
    
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


# ========== eNewick Example Strings ==========


# Simple tree examples
ENEWICK_SIMPLE_TREE = "(A,B);"
ENEWICK_TREE_3_LEAVES = "(A,B,C);"
ENEWICK_TREE_4_LEAVES = "((A,B),(C,D));"
ENEWICK_TREE_WITH_LABELS = "((A,B)int1,(C,D)int2)root;"

# Tree with branch lengths
ENEWICK_TREE_WITH_BRANCH_LENGTHS = "((A:0.5,B:0.3):0.2,(C:0.4,D:0.6):0.1);"

# Hybrid node examples (Extended Newick)
ENEWICK_SINGLE_HYBRID = "((A,B)#H1,#H1);"
ENEWICK_HYBRID_WITH_SIBLING = "((A,B)#H1,(#H1,C));"
ENEWICK_TWO_HYBRIDS = "(((A,B)#H1,C)#H2,(#H1,#H2));"
ENEWICK_HYBRID_COMPLEX = "((((A,B)#H1,C)#H2,D)#H3,(#H1,#H2,#H3));"

# With branch lengths and hybrid (branch length on reference, not definition)
ENEWICK_HYBRID_WITH_BRANCH_LENGTHS = "((A:0.5,B:0.3)#H1,#H1:0.4);"

# With comments (gamma, bootstrap) - comments go on the node before hybrid marker
# Note: Parser doesn't support comments after hybrid markers, so we put them on the node
ENEWICK_WITH_GAMMA = "((A,B)[&gamma=0.6]#H1,#H1);"
ENEWICK_WITH_BOOTSTRAP = "((A,B)[&bootstrap=0.95],(C,D)[&bootstrap=0.87]);"
ENEWICK_WITH_ATTRIBUTES = "((A:0.5,B:0.3)[&bootstrap=0.95][&gamma=0.6]#H1,#H1:0.4);"

# Non-binary (polytomy)
ENEWICK_POLYTOMY = "(A,B,C,D);"
ENEWICK_POLYTOMY_INTERNAL = "((A,B,C)int1,D);"
ENEWICK_HYBRID_THREE_PARENTS = "((A,B)#H1,#H1,#H1);"

# Quoted labels
ENEWICK_QUOTED_LABELS = "('Taxon A','Taxon B','Taxon C');"
ENEWICK_QUOTED_WITH_SPACES = "('Species 1','Species 2','Species 3');"

# Scientific notation
ENEWICK_SCIENTIFIC_NOTATION = "((A:1.5e-2,B:2.3e-1),C:3.4e-2);"

# Larger examples (from literature/common patterns)
ENEWICK_LARGER_TREE = "((((A,B),(C,D)),((E,F),(G,H))),((I,J),(K,L)));"
ENEWICK_LARGER_WITH_HYBRIDS = "((((A,B)#H1,(C,D))#H2,(E,F)),((#H1,#H2),(G,H)));"

# Real-world inspired (5 taxa with hybrid)
ENEWICK_5_TAXA_HYBRID = "((((A,B)#H1,C),D),#H1);"

# Complex nested hybrid
ENEWICK_NESTED_HYBRIDS = "((((A,B)#H1,C)#H2,D),((#H1,E)#H3,#H2,#H3));"


# ========== Helper Functions ==========


def get_all_networks() -> dict[str, DirectedPhyNetwork | SemiDirectedPhyNetwork]:
    """
    Get all network examples as a dictionary.
    
    Returns
    -------
    dict[str, DirectedPhyNetwork | SemiDirectedPhyNetwork]
        Dictionary mapping names to network instances.
    """
    return {
        'empty': empty_network(),
        'single_node': single_node_network(),
        'single_edge': single_edge_network(),
        'two_leaf_star': two_leaf_star(),
        'three_leaf_star': three_leaf_star(),
        'binary_tree_4_leaves': binary_tree_4_leaves(),
        'binary_tree_with_internal': binary_tree_with_internal(),
        'single_hybrid': single_hybrid_network(),
        'single_hybrid_simple': single_hybrid_network_simple(),
        'two_hybrid': two_hybrid_network(),
        'hybrid_parallel_edges': hybrid_with_parallel_edges(),
        'tree_branch_lengths': tree_with_branch_lengths(),
        'network_bootstrap': network_with_bootstrap(),
        'semidirected_simple': semidirected_simple(),
        'semidirected_hybrid': semidirected_with_hybrid(),
        'star_10': star_tree_10_leaves(),
        'star_50': star_tree_50_leaves(),
        'binary_tree_6_levels': balanced_binary_tree(6),
    }


def get_all_enewick_strings() -> dict[str, str]:
    """
    Get all eNewick example strings as a dictionary.
    
    Returns
    -------
    dict[str, str]
        Dictionary mapping names to eNewick strings.
    """
    return {
        'simple_tree': ENEWICK_SIMPLE_TREE,
        'tree_3_leaves': ENEWICK_TREE_3_LEAVES,
        'tree_4_leaves': ENEWICK_TREE_4_LEAVES,
        'tree_with_labels': ENEWICK_TREE_WITH_LABELS,
        'tree_branch_lengths': ENEWICK_TREE_WITH_BRANCH_LENGTHS,
        'single_hybrid': ENEWICK_SINGLE_HYBRID,
        'hybrid_with_sibling': ENEWICK_HYBRID_WITH_SIBLING,
        'two_hybrids': ENEWICK_TWO_HYBRIDS,
        'hybrid_complex': ENEWICK_HYBRID_COMPLEX,
        'hybrid_branch_lengths': ENEWICK_HYBRID_WITH_BRANCH_LENGTHS,
        'with_gamma': ENEWICK_WITH_GAMMA,
        'with_bootstrap': ENEWICK_WITH_BOOTSTRAP,
        'with_attributes': ENEWICK_WITH_ATTRIBUTES,
        'polytomy': ENEWICK_POLYTOMY,
        'polytomy_internal': ENEWICK_POLYTOMY_INTERNAL,
        'hybrid_three_parents': ENEWICK_HYBRID_THREE_PARENTS,
        'quoted_labels': ENEWICK_QUOTED_LABELS,
        'quoted_with_spaces': ENEWICK_QUOTED_WITH_SPACES,
        'scientific_notation': ENEWICK_SCIENTIFIC_NOTATION,
        'larger_tree': ENEWICK_LARGER_TREE,
        'larger_with_hybrids': ENEWICK_LARGER_WITH_HYBRIDS,
        '5_taxa_hybrid': ENEWICK_5_TAXA_HYBRID,
        'nested_hybrids': ENEWICK_NESTED_HYBRIDS,
    }


# ========== Network Generators ==========


def generate_star_tree(n_leaves: int, root_id: int = 1000) -> DirectedPhyNetwork:
    """
    Generate a star tree with n leaves.
    
    Parameters
    ----------
    n_leaves : int
        Number of leaves.
    root_id : int, optional
        ID for root node. By default 1000.
    
    Returns
    -------
    DirectedPhyNetwork
        Star tree network.
    """
    edges = [(root_id, i) for i in range(1, n_leaves + 1)]
    nodes = [(i, {'label': f'Taxon{i}'}) for i in range(1, n_leaves + 1)]
    return DirectedPhyNetwork(edges=edges, nodes=nodes)


def generate_balanced_tree(n_leaves: int, root_id: int = 10000) -> DirectedPhyNetwork:
    """
    Generate a balanced binary tree with approximately n leaves.
    
    Parameters
    ----------
    n_leaves : int
        Target number of leaves (actual may vary slightly).
    root_id : int, optional
        Starting ID for root node. By default 10000.
    
    Returns
    -------
    DirectedPhyNetwork
        Balanced binary tree.
    """
    import math
    depth = math.ceil(math.log2(n_leaves))
    return balanced_binary_tree(depth)

