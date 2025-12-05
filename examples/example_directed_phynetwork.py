"""
Example demonstrating DirectedPhyNetwork usage.

This example shows how to create and work with directed phylogenetic networks,
including trees, hybrid networks, edge attributes, and network properties.

To run this example:
    python examples/example_directed_phynetwork.py

Or from the project root:
    python -m examples.example_directed_phynetwork
"""
import sys
import warnings
from pathlib import Path

# Add src to path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

from phylozoo.core.network.dnetwork.d_phynetwork import DirectedPhyNetwork


def example_basic_tree() -> None:
    """Example 1: Creating a basic phylogenetic tree."""
    print("=" * 60)
    print("Example 1: Basic Phylogenetic Tree")
    print("=" * 60)
    
    # Create a simple tree with 3 leaves
    net = DirectedPhyNetwork(
        edges=[(3, 1), (3, 2), (3, 4)],
        taxa={1: "A", 2: "B", 4: "C"}
    )
    
    print(f"Network: {net}")
    print(f"Root node: {net.root_node}")
    print(f"Leaves: {sorted(net.leaves)}")
    print(f"Taxa: {sorted(net.taxa)}")
    print(f"Is tree: {net.is_tree()}")
    print(f"Number of nodes: {net.number_of_nodes()}")
    print(f"Number of edges: {net.number_of_edges()}")
    print()


def example_tree_with_attributes() -> None:
    """Example 2: Tree with branch lengths and bootstrap values."""
    print("=" * 60)
    print("Example 2: Tree with Edge Attributes")
    print("=" * 60)
    
    # Create a tree with branch lengths and bootstrap support
    net = DirectedPhyNetwork(
        edges=[
            {'u': 5, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
            {'u': 5, 'v': 2, 'branch_length': 0.3, 'bootstrap': 0.87},
            {'u': 5, 'v': 6, 'branch_length': 0.4, 'bootstrap': 0.92},
            {'u': 6, 'v': 3, 'branch_length': 0.2, 'bootstrap': 0.78},
            {'u': 6, 'v': 4, 'branch_length': 0.25, 'bootstrap': 0.81}
        ],
        taxa={1: "A", 2: "B", 3: "C", 4: "D"}
    )
    
    print(f"Network: {net}")
    print(f"Root: {net.root_node}")
    
    # Access edge attributes
    print("\nEdge attributes:")
    for u, v in net.tree_edges:
        bl = net.get_branch_length(u, v)
        bs = net.get_bootstrap(u, v)
        if bl is not None or bs is not None:
            print(f"  Edge ({u}, {v}): branch_length={bl}, bootstrap={bs}")
    print()


def example_hybrid_network() -> None:
    """Example 3: Network with hybrid nodes and gamma values."""
    print("=" * 60)
    print("Example 3: Hybrid Network with Gamma Values")
    print("=" * 60)
    
    # Create a network with a hybrid node
    # Hybrid node 4 receives genetic material from nodes 5 and 6
    net = DirectedPhyNetwork(
        edges=[
            (7, 5), (7, 6),  # Root to tree nodes
            {'u': 5, 'v': 4, 'gamma': 0.6, 'branch_length': 0.3},  # Hybrid edge
            {'u': 6, 'v': 4, 'gamma': 0.4, 'branch_length': 0.2},  # Hybrid edge
            (5, 8), (6, 9),  # Tree nodes also have other children
            (4, 1)  # Hybrid to leaf
        ],
        taxa={1: "A", 8: "B", 9: "C"}
    )
    
    print(f"Network: {net}")
    print(f"Root: {net.root_node}")
    print(f"Hybrid nodes: {net.hybrid_nodes}")
    print(f"Tree nodes: {sorted(net.tree_nodes)}")
    print(f"Is tree: {net.is_tree()}")
    
    # Access gamma values (only on hybrid edges)
    print("\nHybrid edges with gamma values:")
    for u, v in net.hybrid_edges:
        gamma = net.get_gamma(u, v)
        bl = net.get_branch_length(u, v)
        print(f"  Edge ({u}, {v}): gamma={gamma}, branch_length={bl}")
    
    # Verify gamma sum
    hybrid_node = net.hybrid_nodes[0]
    gamma_sum = sum(
        net.get_gamma(u, v) or 0.0
        for u, v in net.hybrid_edges
        if v == hybrid_node
    )
    print(f"\nGamma sum for hybrid node {hybrid_node}: {gamma_sum}")
    print()


def example_complex_network() -> None:
    """Example 4: Complex network with multiple hybrid nodes."""
    print("=" * 60)
    print("Example 4: Complex Network with Multiple Hybrids")
    print("=" * 60)
    
    # Create a more complex network with multiple hybrid nodes
    net = DirectedPhyNetwork(
        edges=[
            # Root structure
            (10, 8), (10, 9),
            # First hybrid node (4)
            (8, 5), (8, 6),
            {'u': 5, 'v': 4, 'gamma': 0.7},
            {'u': 6, 'v': 4, 'gamma': 0.3},
            (5, 1), (6, 2),  # Leaves from first hybrid's parents
            # Second hybrid node (7)
            (9, 11), (9, 12),  # Tree node 9 splits to 11, 12
            {'u': 11, 'v': 7, 'gamma': 0.5},
            {'u': 4, 'v': 7, 'gamma': 0.5},  # Hybrid 4 connects to hybrid 7
            (7, 3),  # Hybrid 7 to leaf
            (11, 13), (11, 14),  # Tree node 11 splits
            (12, 15), (12, 16)  # Tree node 12 splits
        ],
        taxa={1: "A", 2: "B", 3: "C", 13: "D", 14: "E", 15: "F", 16: "G"}
    )
    
    print(f"Network: {net}")
    print(f"Root: {net.root_node}")
    print(f"Hybrid nodes: {sorted(net.hybrid_nodes)}")
    print(f"Tree nodes: {sorted(net.tree_nodes)}")
    print(f"Leaves: {sorted(net.leaves)}")
    print(f"Taxa: {sorted(net.taxa)}")
    
    print("\nHybrid edges:")
    for u, v in sorted(net.hybrid_edges):
        gamma = net.get_gamma(u, v)
        print(f"  ({u}, {v}): gamma={gamma}")
    
    print(f"\nNumber of hybrid edges: {len(net.hybrid_edges)}")
    print(f"Number of tree edges: {len(net.tree_edges)}")
    print()


def example_node_labels() -> None:
    """Example 5: Working with node labels."""
    print("=" * 60)
    print("Example 5: Node Labels")
    print("=" * 60)
    
    # Create a network with labeled internal nodes
    net = DirectedPhyNetwork(
        edges=[(3, 1), (3, 2)],
        taxa={1: "Species_A", 2: "Species_B"},
        internal_node_labels={3: "MRCA"}
    )
    
    print(f"Network: {net}")
    print(f"Taxa: {sorted(net.taxa)}")
    
    # Access labels
    print("\nNode labels:")
    for node in net:
        label = net.get_label(node)
        if label:
            print(f"  Node {node}: '{label}'")
        else:
            print(f"  Node {node}: (no label)")
    
    # Look up nodes by label
    print("\nLookup by label:")
    node_id = net.get_node_id("Species_A")
    print(f"  'Species_A' -> Node {node_id}")
    node_id = net.get_node_id("MRCA")
    print(f"  'MRCA' -> Node {node_id}")
    print()


def example_network_properties() -> None:
    """Example 6: Exploring network topology properties."""
    print("=" * 60)
    print("Example 6: Network Topology Properties")
    print("=" * 60)
    
    net = DirectedPhyNetwork(
        edges=[
            (10, 8), (10, 9),
            (8, 5), (8, 6),
            {'u': 5, 'v': 4, 'gamma': 0.6},
            {'u': 6, 'v': 4, 'gamma': 0.4},
            (5, 1), (6, 2),
            (4, 3)
        ],
        taxa={1: "A", 2: "B", 3: "C", 9: "D"}
    )
    
    print(f"Network: {net}")
    print(f"\nTopology:")
    print(f"  Root node: {net.root_node}")
    print(f"  Leaves: {sorted(net.leaves)}")
    print(f"  Internal nodes: {sorted(net.internal_nodes)}")
    print(f"  Tree nodes: {sorted(net.tree_nodes)}")
    print(f"  Hybrid nodes: {sorted(net.hybrid_nodes)}")
    
    print(f"\nNode degrees:")
    for node in sorted(net):
        in_deg = net.indegree(node)
        out_deg = net.outdegree(node)
        node_type = "root" if node == net.root_node else \
                   "leaf" if node in net.leaves else \
                   "hybrid" if node in net.hybrid_nodes else "tree"
        print(f"  Node {node} ({node_type}): in-degree={in_deg}, out-degree={out_deg}")
    
    print(f"\nIncident edges:")
    for node in sorted(net):
        if node != net.root_node:
            parent_edges = list(net.incident_parent_edges(node))
            print(f"  Node {node} incoming edges: {parent_edges}")
        if node not in net.leaves:
            child_edges = list(net.incident_child_edges(node))
            print(f"  Node {node} outgoing edges: {child_edges}")
    print()


def example_auto_labeling() -> None:
    """Example 7: Automatic labeling of uncovered leaves."""
    print("=" * 60)
    print("Example 7: Automatic Leaf Labeling")
    print("=" * 60)
    
    # Create a network but only label some leaves
    # Uncovered leaves get auto-generated labels
    net = DirectedPhyNetwork(
        edges=[(3, 1), (3, 2), (3, 4), (3, 5)],
        taxa={1: "Species_A", 4: "Species_C"}  # Only label 2 leaves
    )
    
    print(f"Network: {net}")
    print(f"All leaves: {sorted(net.leaves)}")
    print(f"Taxa (all leaf labels): {sorted(net.taxa)}")
    print("\nNote: Leaves 2 and 4 were auto-labeled with their node IDs")
    print()


def example_empty_network() -> None:
    """Example 8: Empty network (with warning)."""
    print("=" * 60)
    print("Example 8: Empty Network")
    print("=" * 60)
    
    # Create an empty network (raises a warning)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        net = DirectedPhyNetwork(edges=[])
        if w:
            print(f"Warning caught: {w[0].message}")
    
    print(f"Empty network: {net}")
    print(f"Number of nodes: {net.number_of_nodes()}")
    print(f"Number of edges: {net.number_of_edges()}")
    print(f"Is valid: {net.validate()}")
    print()


def main() -> None:
    """Run all DirectedPhyNetwork examples."""
    print("\n" + "=" * 60)
    print("DirectedPhyNetwork Examples")
    print("=" * 60 + "\n")
    
    example_basic_tree()
    example_tree_with_attributes()
    example_hybrid_network()
    example_complex_network()
    example_node_labels()
    example_network_properties()
    example_auto_labeling()
    example_empty_network()
    
    print("=" * 60)
    print("All examples completed!")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()

