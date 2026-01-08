"""
Example script demonstrating visualize2 plotting with a larger network.

This script plots a larger network from the test fixtures to better
showcase the combining view layout algorithm.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt

# Add project root to path for test fixtures
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize2 import NetworkStyle, plot_network
from tests.fixtures.directed_networks import (
    LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
    LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB,
    LEVEL_2_DNETWORK_DIAMOND_HYBRID,
)


def main() -> None:
    """Plot larger example networks."""
    # Use a network with multiple hybrids
    networks = [
        ('Two Hybrids Separate', LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE),
        ('Two Hybrids Same Blob', LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB),
        ('Diamond Hybrid', LEVEL_2_DNETWORK_DIAMOND_HYBRID),
    ]

    fig = plt.figure(figsize=(18, 6))
    fig.suptitle('Visualize2: Larger Networks with Combining View', fontsize=16, fontweight='bold')

    for idx, (name, network) in enumerate(networks, 1):
        ax = plt.subplot(1, 3, idx)
        style = NetworkStyle(
            node_color='lightblue',
            leaf_color='lightgreen',
            hybrid_color='salmon',
            edge_color='gray',
            hybrid_edge_color='red',
            edge_width=2.0,
        )
        plot_network(network, style=style, ax=ax, iterations=300)
        ax.set_title(f'{idx}. {name}\n({network.number_of_nodes()} nodes, {network.number_of_edges()} edges)', 
                     fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig('visualize2_large_example.png', dpi=150, bbox_inches='tight')
    print("Plot saved to visualize2_large_example.png")
    plt.show()


if __name__ == '__main__':
    main()

