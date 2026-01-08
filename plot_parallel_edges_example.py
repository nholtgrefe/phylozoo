"""
Example script to demonstrate parallel edge rendering in rectangular layouts.
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from phylozoo.visualize.network_plot import plot_network_with_layout_type
from tests.fixtures.directed_networks import LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES


def main():
    """Plot a network with parallel edges to demonstrate the rendering."""
    # Use a network with many parallel edges
    net = LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES

    print(f"Network: {len(net._graph.nodes)} nodes, {len(net._graph.edges)} edges")
    print(f"Hybrid nodes: {net.hybrid_nodes}")
    print(f"Leaves: {net.leaves}")

    # Count parallel edges
    from phylozoo.visualize.layout.dnetwork import compute_rectangular_dnet_layout

    layout = compute_rectangular_dnet_layout(net)
    parallel_hybrid = sum(
        1
        for route in layout.edge_routes.values()
        if route.edge_type.is_hybrid and route.edge_type.is_parallel
    )
    print(f"Parallel hybrid edges: {parallel_hybrid}")

    # Create plot
    print("\nCreating plot with parallel edge rendering...")
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
    title = f"Network with Parallel Edges\n"
    title += f"Nodes: {len(net._graph.nodes)}, Edges: {len(net._graph.edges)}, "
    title += f"Parallel Hybrid Edges: {parallel_hybrid}"
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
        Line2D([0], [0], color='blue', lw=2, label='Parallel Edge (Downward Curve)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig('parallel_edges_example.png', dpi=150, bbox_inches='tight')
    print("Plot saved to 'parallel_edges_example.png'")
    print("Displaying plot...")
    plt.show()


if __name__ == "__main__":
    main()

