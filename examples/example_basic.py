"""
Basic example demonstrating phylozoo usage.

This example shows how to create and work with basic phylozoo objects.
"""
from phylozoo.core import Split, SplitSystem
from phylozoo.networks import DirectedNetwork
from phylozoo.msa import MSA


def main() -> None:
    """Run basic phylozoo examples."""
    print("=== PhyloZoo Basic Example ===\n")

    # Example 1: Working with splits
    print("1. Creating splits:")
    split1 = Split({1, 2}, {3, 4})
    split2 = Split({1, 3}, {2, 4})
    print(f"   Split 1: {split1}")
    print(f"   Split 2: {split2}")
    print(f"   Split 1 is trivial: {split1.is_trivial()}\n")

    # Example 2: Creating a split system
    print("2. Creating a split system:")
    system = SplitSystem([split1, split2])
    print(f"   Split system with {len(system)} splits\n")

    # Example 3: Working with networks
    print("3. Creating a directed network:")
    network = DirectedNetwork()
    network.add_node("A")
    network.add_node("B")
    network.add_edge("A", "B")
    print(f"   Network: {network}\n")

    # Example 4: Working with MSA
    print("4. Creating an MSA:")
    msa = MSA()
    msa.add_sequence("seq1", "ATCG")
    msa.add_sequence("seq2", "ATCC")
    print(f"   MSA with {len(msa)} sequences\n")

    print("=== Example Complete ===")


if __name__ == "__main__":
    main()

