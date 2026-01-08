"""
Example script to demonstrate the updated rectangular layout computation.
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from phylozoo.visualize.network_plot import plot_network_with_layout_type
from tests.fixtures.directed_networks import (
    LEVEL_2_DNETWORK_NESTED_HYBRIDS,
    LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES,
)


def main():
    """Plot networks with the updated rectangular layout."""
    # Use a network with nested hybrids
    net = LEVEL_2_DNETWORK_NESTED_HYBRIDS

    print(f"Network: {len(net._graph.nodes)} nodes, {len(net._graph.edges)} edges")
    print(f"Hybrid nodes: {net.hybrid_nodes}")
    print(f"Leaves: {net.leaves}")

    # Create plot
    print("\nCreating plot with updated rectangular layout...")
    fig, ax = plt.subplots(figsize=(16, 12))

    plot_network_with_layout_type(
        net,
        layout_type='rectangular',
        ax=ax,
        layer_gap=1.5,
        node_gap=1.2,
        iterations=50,
        seed=42,
        curve_strength=0.4,
    )

    # Add title
    title = f"Updated Rectangular Layout (Dendroscope-style)\n"
    title += f"Nodes: {len(net._graph.nodes)}, Edges: {len(net._graph.edges)}, "
    title += f"Hybrid Nodes: {len(net.hybrid_nodes)}"
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    # Add legend
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            markerfacecolor='lightblue',
            markeredgecolor='darkblue',
            markersize=10,
            label='Tree Node',
        ),
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            markerfacecolor='salmon',
            markeredgecolor='darkred',
            markersize=10,
            label='Hybrid Node',
        ),
        Line2D(
            [0],
            [0],
            marker='s',
            color='w',
            markerfacecolor='lightgreen',
            markeredgecolor='darkgreen',
            markersize=10,
            label='Leaf',
        ),
        Line2D([0], [0], color='gray', lw=2, label='Tree Edge'),
        Line2D([0], [0], color='red', lw=2, label='Hybrid Edge (Curved)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig('updated_layout_example.png', dpi=150, bbox_inches='tight')
    print("Plot saved to 'updated_layout_example.png'")

    # Also create a plot with parallel edges
    print("\nCreating plot with parallel edges...")
    net_parallel = LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES
    fig2, ax2 = plt.subplots(figsize=(16, 12))

    plot_network_with_layout_type(
        net_parallel,
        layout_type='rectangular',
        ax=ax2,
        layer_gap=1.5,
        node_gap=1.2,
        iterations=50,
        seed=42,
        curve_strength=0.4,
    )

    # Count parallel edges
    from phylozoo.visualize.layout.dnetwork import compute_rectangular_dnet_layout

    layout = compute_rectangular_dnet_layout(net_parallel, seed=42)
    parallel_hybrid = sum(
        1
        for route in layout.edge_routes.values()
        if route.edge_type.is_hybrid and route.edge_type.is_parallel
    )

    title2 = f"Updated Layout with Parallel Edges\n"
    title2 += f"Nodes: {len(net_parallel._graph.nodes)}, Edges: {len(net_parallel._graph.edges)}, "
    title2 += f"Parallel Hybrid Edges: {parallel_hybrid}"
    ax2.set_title(title2, fontsize=14, fontweight='bold', pad=20)

    legend_elements2 = [
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            markerfacecolor='lightblue',
            markeredgecolor='darkblue',
            markersize=10,
            label='Tree Node',
        ),
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            markerfacecolor='salmon',
            markeredgecolor='darkred',
            markersize=10,
            label='Hybrid Node',
        ),
        Line2D(
            [0],
            [0],
            marker='s',
            color='w',
            markerfacecolor='lightgreen',
            markeredgecolor='darkgreen',
            markersize=10,
            label='Leaf',
        ),
        Line2D([0], [0], color='gray', lw=2, label='Tree Edge'),
        Line2D([0], [0], color='red', lw=2, label='Hybrid Edge (Curved)'),
        Line2D([0], [0], color='blue', lw=2, label='Parallel Edge (Downward Curve)'),
    ]
    ax2.legend(handles=legend_elements2, loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig('updated_layout_parallel_example.png', dpi=150, bbox_inches='tight')
    print("Plot saved to 'updated_layout_parallel_example.png'")

    print("\nDisplaying plots...")
    plt.show()


if __name__ == "__main__":
    main()

