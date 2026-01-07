"""
Test script to verify all network examples are valid.

This script loads all networks from network_examples.py and validates them.
"""

import sys
from pathlib import Path

# Add examples directory to path
examples_dir = Path(__file__).parent
sys.path.insert(0, str(examples_dir))

from network_examples import get_all_networks, get_all_enewick_strings


def test_all_networks():
    """Test that all network examples are valid."""
    print("Testing all network examples...")
    print("=" * 60)
    
    networks = get_all_networks()
    failed = []
    
    for name, net in networks.items():
        try:
            net.validate()
            print(f"✓ {name:30s} - Valid ({net.number_of_nodes()} nodes, {net.number_of_edges()} edges)")
        except Exception as e:
            print(f"✗ {name:30s} - FAILED: {e}")
            failed.append((name, str(e)))
    
    print("=" * 60)
    if failed:
        print(f"\n{len(failed)} network(s) failed validation:")
        for name, error in failed:
            print(f"  - {name}: {error}")
        return False
    else:
        print(f"\nAll {len(networks)} networks are valid!")
        return True


def test_enewick_strings():
    """Test that all eNewick strings can be parsed."""
    print("\nTesting eNewick string parsing...")
    print("=" * 60)
    
    try:
        from phylozoo.core.network.dnetwork._enewick import parse_enewick
    except ImportError:
        print("Skipping eNewick tests (parser not available)")
        return True
    
    enewick_strings = get_all_enewick_strings()
    failed = []
    
    for name, enewick_str in enewick_strings.items():
        try:
            parsed = parse_enewick(enewick_str)
            print(f"✓ {name:30s} - Parsed ({len(parsed.nodes)} nodes, {len(parsed.edges)} edges)")
        except Exception as e:
            print(f"✗ {name:30s} - FAILED: {e}")
            failed.append((name, str(e)))
    
    print("=" * 60)
    if failed:
        print(f"\n{len(failed)} eNewick string(s) failed parsing:")
        for name, error in failed:
            print(f"  - {name}: {error}")
        return False
    else:
        print(f"\nAll {len(enewick_strings)} eNewick strings parsed successfully!")
        return True


def main():
    """Run all tests."""
    print("Network Examples Validation")
    print("=" * 60)
    print()
    
    networks_ok = test_all_networks()
    enewick_ok = test_enewick_strings()
    
    print()
    print("=" * 60)
    if networks_ok and enewick_ok:
        print("All tests passed!")
    else:
        print("Some tests failed.")
    print("=" * 60)
    
    return networks_ok and enewick_ok


if __name__ == '__main__':
    import sys
    success = main()
    sys.exit(0 if success else 1)

