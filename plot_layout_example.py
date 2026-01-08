#!/usr/bin/env python3
"""
Example script to plot a network using the new RectangularLayout class.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

import matplotlib.pyplot as plt

from phylozoo.visualize.network_plot import plot_network_with_layout_type
from tests.fixtures.directed_networks import LEVEL_2_DNETWORK_NESTED_HYBRIDS




def main():
    """Main function to create and display the plot."""
    # Load network
    net = LEVEL_2_DNETWORK_NESTED_HYBRIDS
    
    print(f"Network: {len(net._graph.nodes)} nodes, {len(net._graph.edges)} edges")
    print(f"Hybrid nodes: {net.hybrid_nodes}")
    print(f"Leaves: {net.leaves}")
    
    # Plot using the new public API
    print("\nCreating plot with rectangular layout...")
    fig, ax = plt.subplots(figsize=(16, 12))
    plot_network_with_layout_type(
        net,
        layout_type='rectangular',
        ax=ax,
        layer_gap=1.5,
        node_gap=1.2,
        iterations=50,
        curve_strength=0.4,
    )
    
    # Add title
    title = f"Level-2 Network with {len(net.hybrid_nodes)} Hybrid Nodes\n"
    title += f"Nodes: {len(net._graph.nodes)}, Edges: {len(net._graph.edges)}"
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue',
               markeredgecolor='darkblue', markersize=10, label='Tree Node'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='salmon',
               markeredgecolor='darkred', markersize=10, label='Hybrid Node'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='lightgreen',
               markeredgecolor='darkgreen', markersize=10, label='Leaf'),
        Line2D([0], [0], color='gray', lw=2, label='Tree Edge'),
        Line2D([0], [0], color='red', lw=2, label='Hybrid Edge'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('layout_example.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved to 'layout_example.png'")
    print("To display the plot, uncomment plt.show() or run in an interactive environment")
    # plt.show()  # Uncomment to display interactively


if __name__ == '__main__':
    main()

