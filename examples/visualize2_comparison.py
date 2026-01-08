"""
Comparison script: old visualize vs new visualize2.

This script plots the same network using both implementations to compare results.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize2 import NetworkStyle, plot_network
from phylozoo.visualize.layout.dnetwork import compute_combining_dnet_layout, render_combining_dnet_layout
from phylozoo.visualize.network_plot import _style_combining_layout
from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID


def main() -> None:
    """Compare old and new implementations."""
    network = LEVEL_1_DNETWORK_SINGLE_HYBRID

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle('Comparison: Old visualize vs New visualize2', fontsize=16, fontweight='bold')

    # Old implementation
    ax1.set_title('Old visualize module', fontsize=12, fontweight='bold')
    layout_old = compute_combining_dnet_layout(network, style='cladogram', iterations=200)
    render_elements = render_combining_dnet_layout(layout_old, ax1)
    _style_combining_layout(
        layout_old, ax1, render_elements,
        node_color='lightblue',
        leaf_color='lightgreen',
        hybrid_color='salmon',
        edge_color='gray',
        hybrid_edge_color='red',
        edge_width=2.0,
        with_labels=True,
    )
    ax1.set_aspect('equal')
    ax1.axis('off')

    # New implementation
    ax2.set_title('New visualize2 module', fontsize=12, fontweight='bold')
    style = NetworkStyle(
        node_color='lightblue',
        leaf_color='lightgreen',
        hybrid_color='salmon',
        edge_color='gray',
        hybrid_edge_color='red',
        edge_width=2.0,
    )
    plot_network(network, style=style, ax=ax2, iterations=200)

    plt.tight_layout()
    plt.savefig('visualize2_comparison.png', dpi=150, bbox_inches='tight')
    print("Comparison plot saved to visualize2_comparison.png")
    plt.show()


if __name__ == '__main__':
    main()

