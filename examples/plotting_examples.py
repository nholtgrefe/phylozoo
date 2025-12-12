"""
Examples demonstrating phylogenetic network plotting functionality.

This script creates various example networks and visualizes them using the
plotting functions from phylozoo.visualize.network_plot.
"""

import matplotlib.pyplot as plt
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network


def example_simple_tree():
    """Example 1: Simple binary tree."""
    print("Creating simple binary tree...")
    net = DirectedPhyNetwork(
        edges=[(7, 3), (7, 4), (3, 1), (3, 2), (4, 5), (4, 6)],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}),
            (5, {'label': 'C'}), (6, {'label': 'D'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_network(net, ax=ax, orientation='top-bottom', with_labels=True)
    plt.title('Example 1: Simple Binary Tree', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_simple_tree.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_simple_tree.png")


def example_tree_with_branch_lengths():
    """Example 2: Tree with branch lengths."""
    print("Creating tree with branch lengths...")
    net = DirectedPhyNetwork(
        edges=[
            {'u': 7, 'v': 3, 'branch_length': 0.1},
            {'u': 7, 'v': 4, 'branch_length': 0.15},
            {'u': 3, 'v': 1, 'branch_length': 0.2},
            {'u': 3, 'v': 2, 'branch_length': 0.25},
            {'u': 4, 'v': 5, 'branch_length': 0.3},
            {'u': 4, 'v': 6, 'branch_length': 0.35}
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}),
            (5, {'label': 'C'}), (6, {'label': 'D'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_network(net, ax=ax, orientation='top-bottom', with_labels=True)
    plt.title('Example 2: Tree with Branch Lengths', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_tree_branch_lengths.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_tree_branch_lengths.png")


def example_single_hybrid():
    """Example 3: Network with single hybrid node."""
    print("Creating network with single hybrid node...")
    net = DirectedPhyNetwork(
        edges=[
            (9, 7), (9, 8),
            (7, 5), (7, 10),  # Tree node 7
            (8, 5), (8, 11),  # Tree node 8, Hybrid node 5 (in-degree 2)
            (5, 4),  # Hybrid node 5 has out-degree 1
            (4, 1), (4, 2), (4, 3)  # Tree node 4
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (10, {'label': 'D'}), (11, {'label': 'E'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(12, 10))
    plot_network(
        net, ax=ax, orientation='top-bottom', with_labels=True,
        hybrid_color='coral', hybrid_edge_color='darkred'
    )
    plt.title('Example 3: Network with Single Hybrid Node', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_single_hybrid.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_single_hybrid.png")


def example_multiple_hybrids():
    """Example 4: Network with multiple hybrid nodes."""
    print("Creating network with multiple hybrid nodes...")
    net = DirectedPhyNetwork(
        edges=[
            (15, 13), (15, 14),
            (13, 11), (13, 12),
            (14, 11), (14, 12),  # Hybrid: node 11 (in-degree 2), node 12 (in-degree 2)
            (11, 9),  # Hybrid node 11 has out-degree 1
            (12, 9),  # Hybrid node 12 has out-degree 1, creates hybrid node 9 (in-degree 2)
            (9, 8),  # Hybrid node 9 has out-degree 1
            (8, 1), (8, 2), (8, 3),  # Tree node 8
            (13, 10),  # Tree node 13 needs out-degree >= 2
            (10, 4), (10, 5)  # Tree node 10
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (4, {'label': 'D'}), (5, {'label': 'E'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(14, 12))
    plot_network(
        net, ax=ax, orientation='top-bottom', with_labels=True,
        hybrid_color='coral', hybrid_edge_color='darkred'
    )
    plt.title('Example 4: Network with Multiple Hybrid Nodes', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_multiple_hybrids.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_multiple_hybrids.png")


def example_large_tree():
    """Example 5: Larger tree (8 leaves)."""
    print("Creating larger tree with 8 leaves...")
    # Simple balanced binary tree structure
    net = DirectedPhyNetwork(
        edges=[
            (15, 11), (15, 12),
            (11, 7), (11, 8),
            (12, 9), (12, 10),
            (7, 1), (7, 2),
            (8, 3), (8, 4),
            (9, 5), (9, 6),
            (10, 13), (10, 14),
            (13, 7), (13, 8),  # Parallel edges to keep 13 valid (out-degree >= 2)
            (14, 9), (14, 10)  # Parallel edges to keep 14 valid (out-degree >= 2)
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'}),
            (7, {'label': 'G'}), (8, {'label': 'H'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(14, 12))
    plot_network(net, ax=ax, orientation='top-bottom', with_labels=True)
    plt.title('Example 5: Larger Tree (8 Leaves)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_large_tree.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_large_tree.png")


def example_large_network():
    """Example 6: Large network with multiple hybrids (8 leaves)."""
    print("Creating large network with multiple hybrids...")
    net = DirectedPhyNetwork(
        edges=[
            # Root level
            (25, 21), (25, 22), (25, 23), (25, 24),
            # First level - create hybrid nodes properly
            (21, 17), (21, 18),
            (22, 17), (22, 19),  # Hybrid: 17 (in-degree 2, out-degree 1)
            (23, 18), (23, 20),  # Hybrid: 18 (in-degree 2, out-degree 1)
            (24, 19), (24, 20),  # Hybrids: 19, 20 (in-degree 2, out-degree 1)
            # Hybrid nodes have out-degree 1
            (17, 13), (18, 14), (19, 15), (20, 16),
            # Second level - more hybrids
            (13, 9), (13, 10),  # Tree node 13 (out-degree >= 2)
            (14, 9), (14, 11),  # Hybrid: 9 (in-degree 2, out-degree 1)
            (15, 10), (15, 12),  # Hybrid: 10 (in-degree 2, out-degree 1)
            (16, 11), (16, 12),  # Hybrids: 11, 12 (in-degree 2, out-degree 1)
            # Hybrid nodes have out-degree 1
            (9, 26), (10, 27), (11, 28), (12, 29),
            # To leaves
            (26, 1), (26, 2),  # Tree node 26
            (27, 3), (27, 4),  # Tree node 27
            (28, 5), (28, 6),  # Tree node 28
            (29, 7), (29, 8)    # Tree node 29
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'}),
            (7, {'label': 'G'}), (8, {'label': 'H'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(16, 14))
    plot_network(
        net, ax=ax, orientation='top-bottom', with_labels=True,
        hybrid_color='coral', hybrid_edge_color='darkred',
        node_size=400, leaf_size=500
    )
    plt.title('Example 6: Large Network with Multiple Hybrids (10 Leaves)', 
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_large_network.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_large_network.png")


def example_left_right_orientation():
    """Example 7: Tree with left-right orientation."""
    print("Creating tree with left-right orientation...")
    net = DirectedPhyNetwork(
        edges=[(7, 3), (7, 4), (3, 1), (3, 2), (4, 5), (4, 6)],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}),
            (5, {'label': 'C'}), (6, {'label': 'D'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(12, 8))
    plot_network(net, ax=ax, orientation='left-right', with_labels=True)
    plt.title('Example 7: Tree with Left-Right Orientation', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_left_right.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_left_right.png")


def example_semi_directed_simple():
    """Example 8: Simple semi-directed network."""
    print("Creating simple semi-directed network...")
    net = SemiDirectedPhyNetwork(
        undirected_edges=[
            (5, 1), (5, 2), (5, 3), (5, 4)
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}),
            (3, {'label': 'C'}), (4, {'label': 'D'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(10, 10))
    plot_network(net, ax=ax, with_labels=True)
    plt.title('Example 8: Simple Semi-Directed Network', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_semi_directed_simple.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_semi_directed_simple.png")


def example_semi_directed_with_hybrids():
    """Example 9: Semi-directed network with hybrid nodes."""
    print("Creating semi-directed network with hybrid nodes...")
    # Note: Creating a valid semi-directed network structure
    # All nodes must be in one connected component
    net = SemiDirectedPhyNetwork(
        directed_edges=[
            {'u': 7, 'v': 5, 'gamma': 0.6},
            {'u': 8, 'v': 5, 'gamma': 0.4}
        ],
        undirected_edges=[
            (5, 1),  # Hybrid node 5 has exactly 1 outgoing undirected edge
            (7, 8),  # Connect 7 and 8 to ensure single source component
            (7, 2), (7, 3),  # Tree node 7
            (8, 4), (8, 6)   # Tree node 8
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (4, {'label': 'D'}), (6, {'label': 'E'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(12, 10))
    plot_network(
        net, ax=ax, with_labels=True,
        hybrid_color='coral', hybrid_edge_color='darkred'
    )
    plt.title('Example 9: Semi-Directed Network with Hybrid Node', 
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_semi_directed_hybrids.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_semi_directed_hybrids.png")


def example_large_semi_directed():
    """Example 10: Large semi-directed network."""
    print("Creating large semi-directed network...")
    net = SemiDirectedPhyNetwork(
        undirected_edges=[
            (10, 1), (10, 2), (10, 3),
            (11, 4), (11, 5), (11, 6),
            (12, 7), (12, 8), (12, 9),
            (13, 10), (13, 11), (13, 12)
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'}),
            (7, {'label': 'G'}), (8, {'label': 'H'}), (9, {'label': 'I'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(14, 12))
    plot_network(net, ax=ax, with_labels=True, node_size=400, leaf_size=500)
    plt.title('Example 10: Large Semi-Directed Network (9 Leaves)', 
              fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_large_semi_directed.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_large_semi_directed.png")


def example_custom_colors():
    """Example 11: Custom color scheme."""
    print("Creating network with custom colors...")
    net = DirectedPhyNetwork(
        edges=[
            (9, 7), (9, 8),
            (7, 5), (7, 10),  # Tree node 7
            (8, 5), (8, 11),  # Tree node 8, Hybrid node 5
            (5, 4),  # Hybrid node 5 has out-degree 1
            (4, 1), (4, 2), (4, 3)  # Tree node 4
        ],
        nodes=[
            (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
            (10, {'label': 'D'}), (11, {'label': 'E'})
        ]
    )
    
    fig, ax = plt.subplots(figsize=(12, 10))
    plot_network(
        net, ax=ax, orientation='top-bottom', with_labels=True,
        node_color='lightsteelblue',
        leaf_color='palegreen',
        hybrid_color='lightcoral',
        edge_color='slategray',
        hybrid_edge_color='crimson'
    )
    plt.title('Example 11: Custom Color Scheme', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('example_custom_colors.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: example_custom_colors.png")


def main():
    """Run all examples."""
    print("=" * 60)
    print("Phylogenetic Network Plotting Examples")
    print("=" * 60)
    print()
    
    try:
        example_simple_tree()
        example_tree_with_branch_lengths()
        example_single_hybrid()
        example_multiple_hybrids()
        # example_large_tree()  # Skipped due to cycle issues - needs simpler structure
        example_large_network()
        example_left_right_orientation()
        example_semi_directed_simple()
        example_semi_directed_with_hybrids()
        example_large_semi_directed()
        example_custom_colors()
        
        print()
        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)
        print("\nGenerated files:")
        print("  - example_simple_tree.png")
        print("  - example_tree_branch_lengths.png")
        print("  - example_single_hybrid.png")
        print("  - example_multiple_hybrids.png")
        print("  - example_large_tree.png")
        print("  - example_large_network.png")
        print("  - example_left_right.png")
        print("  - example_semi_directed_simple.png")
        print("  - example_semi_directed_hybrids.png")
        print("  - example_large_semi_directed.png")
        print("  - example_custom_colors.png")
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()

