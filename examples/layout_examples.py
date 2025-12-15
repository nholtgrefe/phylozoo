"""
Examples demonstrating the new optimal layout functionality.

This script shows how to use the new layout-based plotting functions
for DirectedPhyNetwork and SemiDirectedPhyNetwork with crossing
minimization and optimal positioning.
"""

import matplotlib.pyplot as plt

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network_with_layout

# Import reusable network structures
from .network_examples import (
    binary_tree_4_leaves,
    single_hybrid_network,
    tree_with_branch_lengths,
    semidirected_with_hybrid,
    two_hybrid_network,
)


def example_simple_tree():
    """Example: Simple binary tree with optimal layout."""
    print("Creating simple binary tree example...")
    
    net = binary_tree_4_leaves()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_network_with_layout(
        net,
        ax=ax,
        orientation='top-bottom',
        use_branch_lengths=False,
    )
    plt.tight_layout()
    plt.savefig('example_simple_tree_layout.png', dpi=150, bbox_inches='tight')
    print("  Saved: example_simple_tree_layout.png")
    plt.close()


def example_tree_with_branch_lengths():
    """Example: Tree with branch lengths enabled."""
    print("Creating tree with branch lengths example...")
    
    net = tree_with_branch_lengths()
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_network_with_layout(
        net,
        ax=ax,
        orientation='top-bottom',
        use_branch_lengths=True,  # Enable branch length scaling
    )
    plt.tight_layout()
    plt.savefig('example_tree_with_branch_lengths.png', dpi=150, bbox_inches='tight')
    print("  Saved: example_tree_with_branch_lengths.png")
    plt.close()


def example_network_with_hybrid():
    """Example: Network with hybrid node and crossing minimization."""
    print("Creating network with hybrid node example...")
    
    net = single_hybrid_network()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    plot_network_with_layout(
        net,
        ax=ax,
        orientation='top-bottom',
        hybrid_color='coral',
        hybrid_edge_color='darkred',
    )
    plt.tight_layout()
    plt.savefig('example_network_with_hybrid.png', dpi=150, bbox_inches='tight')
    print("  Saved: example_network_with_hybrid.png")
    plt.close()


def example_left_right_orientation():
    """Example: Left-right orientation for directed network."""
    print("Creating left-right orientation example...")
    
    net = DirectedPhyNetwork(
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
    
    fig, ax = plt.subplots(figsize=(12, 8))
    plot_network_with_layout(
        net,
        ax=ax,
        orientation='left-right',  # Root on left, leaves on right
    )
    plt.tight_layout()
    plt.savefig('example_left_right_orientation.png', dpi=150, bbox_inches='tight')
    print("  Saved: example_left_right_orientation.png")
    plt.close()


def example_semidirected_network():
    """Example: Semi-directed network with leaf repulsion."""
    print("Creating semi-directed network example...")
    
    net = semidirected_with_hybrid()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    plot_network_with_layout(
        net,
        ax=ax,
        leaf_color='lightgreen',
        hybrid_color='salmon',
    )
    plt.tight_layout()
    plt.savefig('example_semidirected_network.png', dpi=150, bbox_inches='tight')
    print("  Saved: example_semidirected_network.png")
    plt.close()


def example_complex_network():
    """Example: More complex network with multiple hybrids."""
    print("Creating complex network example...")
    
    net = two_hybrid_network()
    
    fig, ax = plt.subplots(figsize=(14, 12))
    plot_network_with_layout(
        net,
        ax=ax,
        orientation='top-bottom',
        node_size=400,
        leaf_size=500,
        edge_width=1.5,
    )
    plt.tight_layout()
    plt.savefig('example_complex_network.png', dpi=150, bbox_inches='tight')
    print("  Saved: example_complex_network.png")
    plt.close()


def main():
    """Run all examples."""
    print("=" * 60)
    print("Layout Examples - Optimal Network Plotting")
    print("=" * 60)
    print()
    
    try:
        example_simple_tree()
        example_tree_with_branch_lengths()
        example_network_with_hybrid()
        example_left_right_orientation()
        example_semidirected_network()
        example_complex_network()
        
        print()
        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()

