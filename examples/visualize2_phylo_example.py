"""
Example script demonstrating the phylo layout in visualize2.

This script plots networks using the new phylo layout algorithm.
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
    LEVEL_1_DNETWORK_SINGLE_HYBRID,
    LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
    LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB,
)


def main() -> None:
    """Plot example networks with phylo layout."""
    networks = [
        ('Single Hybrid', LEVEL_1_DNETWORK_SINGLE_HYBRID),
        ('Two Hybrids Separate', LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE),
        ('Two Hybrids Same Blob', LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB),
    ]

    fig = plt.figure(figsize=(18, 12))
    fig.suptitle('Visualize2: Phylo Layout Examples', fontsize=16, fontweight='bold')

    style = NetworkStyle(
        node_color='lightblue',
        leaf_color='lightgreen',
        hybrid_color='salmon',
        edge_color='gray',
        hybrid_edge_color='red',
        edge_width=2.0,
    )

    # Top row: TD (top-down) direction
    for idx, (name, network) in enumerate(networks, 1):
        ax = plt.subplot(2, 3, idx)
        plot_network(
            network,
            layout='phylo',
            style=style,
            ax=ax,
            direction='TD',
            trials=1000,
        )
        ax.set_title(
            f'{idx}. {name} (TD)\n({network.number_of_nodes()} nodes, {network.number_of_edges()} edges)',
            fontsize=12,
            fontweight='bold',
        )

    # Bottom row: LR (left-right) direction
    for idx, (name, network) in enumerate(networks, 4):
        ax = plt.subplot(2, 3, idx)
        plot_network(
            network,
            layout='phylo',
            style=style,
            ax=ax,
            direction='LR',
            trials=1000,
        )
        ax.set_title(
            f'{idx-3}. {name} (LR)',
            fontsize=12,
            fontweight='bold',
        )

    plt.tight_layout()
    plt.savefig('visualize2_phylo_example.png', dpi=150, bbox_inches='tight')
    print("Plot saved to visualize2_phylo_example.png")
    plt.show()


if __name__ == '__main__':
    main()

