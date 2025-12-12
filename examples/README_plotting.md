# Phylogenetic Network Plotting Examples

This directory contains example scripts demonstrating the plotting functionality for phylogenetic networks.

## Running the Examples

To generate all example plots, run:

```bash
python plotting_examples.py
```

This will create 11 PNG files demonstrating various network types and plotting features.

## Examples Included

1. **example_simple_tree.png** - Simple binary tree (4 leaves)
2. **example_tree_branch_lengths.png** - Tree with branch length annotations
3. **example_single_hybrid.png** - Network with a single hybrid node (reticulation)
4. **example_multiple_hybrids.png** - Network with multiple hybrid nodes
5. **example_large_tree.png** - Larger tree structure (8 leaves)
6. **example_large_network.png** - Large network with multiple hybrids (8 leaves, 8 hybrid nodes)
7. **example_left_right.png** - Tree with left-to-right orientation
8. **example_semi_directed_simple.png** - Simple semi-directed network
9. **example_semi_directed_hybrids.png** - Semi-directed network with hybrid nodes
10. **example_large_semi_directed.png** - Large semi-directed network (9 leaves)
11. **example_custom_colors.png** - Network with custom color scheme

## Features Demonstrated

- **Hierarchical layouts** for directed networks (top-bottom and left-right)
- **Planar layouts** for semi-directed networks with labels on outside
- **Hybrid node visualization** (highlighted in different colors)
- **Leaf node highlighting** (green by default)
- **Parallel edge handling** (curved arcs)
- **Customizable colors and sizes**
- **Label positioning** (automatic for directed, outside for semi-directed)

## Usage

```python
from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.visualize.network_plot import plot_network
import matplotlib.pyplot as plt

# Create a network
net = DirectedPhyNetwork(
    edges=[(3, 1), (3, 2)],
    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
)

# Plot it
fig, ax = plt.subplots(figsize=(10, 8))
plot_network(net, ax=ax, with_labels=True)
plt.savefig('my_network.png', dpi=150, bbox_inches='tight')
plt.close()
```

## Customization Options

The `plot_network()` function supports many customization options:

- `orientation`: 'top-bottom' or 'left-right' (for directed networks)
- `node_color`: Color for internal nodes
- `leaf_color`: Color for leaf nodes
- `hybrid_color`: Color for hybrid nodes
- `edge_color`: Color for regular edges
- `hybrid_edge_color`: Color for hybrid edges
- `node_size`: Size of internal nodes
- `leaf_size`: Size of leaf nodes
- `edge_width`: Width of edges
- `with_labels`: Whether to show node labels
- `label_offset`: Offset for labels from nodes

See the examples for more details!

