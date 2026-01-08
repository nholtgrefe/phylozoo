"""
Example script demonstrating visualize2 plotting with different styling options.

This script plots a network from the test fixtures with various styling
configurations to showcase the visualize2 module capabilities.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt

# Add project root to path for test fixtures
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize2 import NetworkStyle, plot_network
from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID


def main() -> None:
    """Plot example network with different styling options."""
    # Use a network with hybrid nodes for more interesting visualization
    network = LEVEL_1_DNETWORK_SINGLE_HYBRID

    # Create figure with subplots for different styling options
    fig = plt.figure(figsize=(16, 12))
    fig.suptitle('Visualize2 Example: Different Styling Options', fontsize=16, fontweight='bold')

    # 1. Default style
    ax1 = plt.subplot(2, 3, 1)
    plot_network(network, ax=ax1)
    ax1.set_title('1. Default Style', fontsize=12, fontweight='bold')

    # 2. Custom colors
    ax2 = plt.subplot(2, 3, 2)
    custom_style = NetworkStyle(
        node_color='steelblue',
        leaf_color='forestgreen',
        hybrid_color='coral',
        edge_color='slategray',
        hybrid_edge_color='crimson',
    )
    plot_network(network, style=custom_style, ax=ax2)
    ax2.set_title('2. Custom Colors', fontsize=12, fontweight='bold')

    # 3. Larger nodes and thicker edges
    ax3 = plt.subplot(2, 3, 3)
    large_style = NetworkStyle(
        node_size=800.0,
        leaf_size=1000.0,
        edge_width=3.0,
    )
    plot_network(network, style=large_style, ax=ax3)
    ax3.set_title('3. Larger Nodes & Thicker Edges', fontsize=12, fontweight='bold')

    # 4. No labels
    ax4 = plt.subplot(2, 3, 4)
    no_labels_style = NetworkStyle(with_labels=False)
    plot_network(network, style=no_labels_style, ax=ax4)
    ax4.set_title('4. No Labels', fontsize=12, fontweight='bold')

    # 5. Phylogram layout (with branch lengths if available)
    ax5 = plt.subplot(2, 3, 5)
    phylogram_style = NetworkStyle(
        node_color='lightblue',
        leaf_color='lightgreen',
        hybrid_color='salmon',
    )
    plot_network(
        network,
        style=phylogram_style,
        ax=ax5,
        layout='combining',  # Layout type
        use_branch_lengths=False,  # Network may not have branch lengths
    )
    ax5.set_title('5. Phylogram Layout', fontsize=12, fontweight='bold')

    # 6. Minimalist style
    ax6 = plt.subplot(2, 3, 6)
    minimalist_style = NetworkStyle(
        node_color='white',
        leaf_color='white',
        hybrid_color='lightgray',
        edge_color='black',
        hybrid_edge_color='darkgray',
        edge_width=1.5,
        node_size=300.0,
        leaf_size=300.0,
        with_labels=True,
        label_font_size=8.0,
    )
    plot_network(network, style=minimalist_style, ax=ax6)
    ax6.set_title('6. Minimalist Style', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig('visualize2_example.png', dpi=150, bbox_inches='tight')
    print("Plot saved to visualize2_example.png")
    plt.show()


if __name__ == '__main__':
    main()

