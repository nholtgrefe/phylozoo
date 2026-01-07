# Network Examples Module

This module (`network_examples.py`) provides a centralized collection of reusable network structures for testing, examples, and demonstrations.

## Purpose

Instead of recreating network structures in each test or example file, you can import pre-defined networks from this module. This ensures:
- **Consistency**: Same network structures used across tests
- **Maintainability**: Update network definitions in one place
- **Reusability**: Easy to find and use common network patterns
- **Validation**: All networks are validated and tested

## Usage

### Importing Networks

```python
from examples.network_examples import (
    binary_tree_4_leaves,
    single_hybrid_network,
    tree_with_branch_lengths,
    semidirected_with_hybrid,
)

# Use in your code
net = binary_tree_4_leaves()
plot_network(net)
```

### Getting All Networks

```python
from examples.network_examples import get_all_networks

networks = get_all_networks()
for name, net in networks.items():
    print(f"{name}: {net.number_of_nodes()} nodes")
```

### Using eNewick Strings

```python
from examples.network_examples import get_all_enewick_strings
from phylozoo.core.network.dnetwork._enewick import parse_enewick

enewick_strings = get_all_enewick_strings()
for name, enewick_str in enewick_strings.items():
    parsed = parse_enewick(enewick_str)
    print(f"{name}: {len(parsed.nodes)} nodes")
```

## Available Networks

### Minimal Networks
- `empty_network()` - Empty network (no nodes, no edges)
- `single_node_network()` - Single-node network
- `single_edge_network()` - Network with single edge (root -> leaf)

### Simple Trees
- `two_leaf_star()` - Star tree with root and two leaves
- `three_leaf_star()` - Star tree with root and three leaves
- `binary_tree_4_leaves()` - Balanced binary tree with 4 leaves
- `binary_tree_with_internal()` - Binary tree with internal node

### Networks with Hybrid Nodes
- `single_hybrid_network()` - Network with single hybrid node
- `single_hybrid_network_simple()` - Simpler network with single hybrid
- `two_hybrid_network()` - Network with two hybrid nodes
- `hybrid_with_parallel_edges()` - Network with hybrid and parallel edges

### Networks with Attributes
- `tree_with_branch_lengths()` - Tree with branch lengths on all edges
- `network_with_bootstrap()` - Network with bootstrap support values

### Semi-Directed Networks
- `semidirected_simple()` - Simple semi-directed network (tree-like)
- `semidirected_with_hybrid()` - Semi-directed network with hybrid node

### Large Networks
- `star_tree_10_leaves()` - Star tree with 10 leaves
- `star_tree_50_leaves()` - Star tree with 50 leaves
- `balanced_binary_tree(depth=6)` - Balanced binary tree generator

### Network Generators
- `generate_star_tree(n_leaves, root_id=1000)` - Generate star tree with n leaves
- `generate_balanced_tree(n_leaves, root_id=10000)` - Generate balanced binary tree

## eNewick Examples

The module also includes a collection of eNewick (Extended Newick) format strings:

- Simple trees: `ENEWICK_SIMPLE_TREE`, `ENEWICK_TREE_3_LEAVES`, etc.
- Trees with branch lengths: `ENEWICK_TREE_WITH_BRANCH_LENGTHS`
- Hybrid networks: `ENEWICK_SINGLE_HYBRID`, `ENEWICK_TWO_HYBRIDS`, etc.
- With attributes: `ENEWICK_WITH_GAMMA`, `ENEWICK_WITH_BOOTSTRAP`
- Larger examples: `ENEWICK_LARGER_TREE`, `ENEWICK_LARGER_WITH_HYBRIDS`

Use `get_all_enewick_strings()` to get all eNewick examples as a dictionary.

## Testing

Run the validation script to ensure all networks are valid:

```bash
python examples/test_network_examples.py
```

This will:
1. Validate all network structures
2. Test parsing of all eNewick strings
3. Report any failures

## Adding New Networks

When adding new network structures:

1. **Follow naming conventions**: Use descriptive names like `network_with_feature()`
2. **Add docstrings**: Document what the network represents
3. **Validate**: Ensure the network passes validation
4. **Add to get_all_networks()**: Include in the dictionary for easy access
5. **Test**: Run `test_network_examples.py` to verify

## Notes

- All networks are validated according to `DirectedPhyNetwork`, `MixedPhyNetwork`, and `SemiDirectedPhyNetwork` rules
- Networks use consistent node ID schemes (avoid conflicts when combining)
- Leaf nodes are always labeled
- Internal nodes may or may not have labels

