"""
Test fixtures for phylogenetic networks.

This package provides a comprehensive collection of pre-built networks
for use in tests. Networks are organized by type and include metadata
documenting their properties.

Usage
-----
Import networks directly from the module:
    
    from tests.fixtures.directed_networks import TREE_SMALL_BINARY
    from tests.fixtures.sd_networks import SINGLE_HYBRID_LEVEL_1

Or access all networks in a category:
    
    from tests.fixtures import directed_networks
    trees = directed_networks.get_networks_by_category('trees')

Or filter by properties:
    
    level_2_nets = directed_networks.get_networks_with_property(level=2)
"""

# Re-export for convenience
from . import directed_networks
from . import sd_networks

__all__ = ['directed_networks', 'sd_networks']


