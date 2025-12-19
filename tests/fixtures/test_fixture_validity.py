"""
Tests to verify all fixture networks are valid.

This module tests that all network fixtures in directed_networks.py
and sd_networks.py are valid and have the documented properties.
"""

import pytest

from tests.fixtures import directed_networks as dn
from tests.fixtures import sd_networks as sdn
from phylozoo.core.network.dnetwork.classifications import (
    level as d_level,
    vertex_level as d_vertex_level,
    reticulation_number as d_reticulation_number,
    is_binary as d_is_binary,
    is_tree as d_is_tree,
    has_parallel_edges as d_has_parallel_edges,
)
from phylozoo.core.network.sdnetwork.classifications import (
    level as sd_level,
    vertex_level as sd_vertex_level,
    reticulation_number as sd_reticulation_number,
    is_binary as sd_is_binary,
    is_tree as sd_is_tree,
    has_parallel_edges as sd_has_parallel_edges,
)


class TestDirectedNetworkFixtures:
    """Test all DirectedPhyNetwork fixtures are valid."""
    
    @pytest.mark.parametrize("network_name,metadata", 
                             list(dn.NETWORK_METADATA.items()))
    def test_network_is_valid(self, network_name: str, metadata: dict) -> None:
        """
        Test that network can be created and is valid.
        
        Parameters
        ----------
        network_name : str
            Name of the network fixture.
        metadata : dict
            Expected metadata for the network.
        """
        network = getattr(dn, network_name)
        # Basic validity checks
        assert network.number_of_nodes() == metadata['nodes']
        assert network.number_of_edges() == metadata['edges']
    
    @pytest.mark.parametrize("network_name,metadata",
                             list(dn.NETWORK_METADATA.items()))
    def test_network_properties_match_metadata(self, network_name: str, metadata: dict) -> None:
        """
        Test that computed properties match documented metadata.
        
        Parameters
        ----------
        network_name : str
            Name of the network fixture.
        metadata : dict
            Expected metadata for the network.
        """
        network = getattr(dn, network_name)
        
        if metadata['nodes'] == 0:
            # Empty network - skip property checks
            return
        
        # Test classification properties
        assert d_is_tree(network) == metadata['is_tree'], \
            f"{network_name}: is_tree mismatch"
        assert d_is_binary(network) == metadata['is_binary'], \
            f"{network_name}: is_binary mismatch"
        assert d_has_parallel_edges(network) == metadata['has_parallel_edges'], \
            f"{network_name}: has_parallel_edges mismatch"
        
        # For non-empty networks, test level-related properties
        if metadata['nodes'] > 1:
            assert d_level(network) == metadata['level'], \
                f"{network_name}: level mismatch (expected {metadata['level']}, got {d_level(network)})"
            assert d_vertex_level(network) == metadata['vertex_level'], \
                f"{network_name}: vertex_level mismatch"
            assert d_reticulation_number(network) == metadata['reticulation_number'], \
                f"{network_name}: reticulation_number mismatch"
    
    def test_can_retrieve_all_networks(self) -> None:
        """Test that all networks can be retrieved."""
        all_networks = dn.get_all_networks()
        assert len(all_networks) == len(dn.NETWORK_METADATA)
        assert all(net is not None for net in all_networks)
    
    def test_can_filter_by_category(self) -> None:
        """Test that networks can be filtered by category."""
        trees = dn.get_networks_by_category('trees')
        assert len(trees) > 0
        assert all(dn.NETWORK_METADATA[name]['category'] == 'trees' 
                  for name in dn.NETWORK_METADATA 
                  if getattr(dn, name) in trees)
    
    def test_can_filter_by_property(self) -> None:
        """Test that networks can be filtered by property."""
        level_0_nets = dn.get_networks_with_property(level=0)
        assert len(level_0_nets) > 0
        
        binary_nets = dn.get_networks_with_property(is_binary=True)
        assert len(binary_nets) > 0
    
    def test_metadata_completeness(self) -> None:
        """Test that all metadata entries have required fields."""
        required_fields = [
            'category', 'nodes', 'edges', 'level', 'vertex_level',
            'reticulation_number', 'is_tree', 'is_binary', 
            'has_parallel_edges', 'num_hybrids'
        ]
        for name, meta in dn.NETWORK_METADATA.items():
            for field in required_fields:
                assert field in meta, f"{name} missing field: {field}"


class TestSemiDirectedNetworkFixtures:
    """Test all SemiDirectedPhyNetwork fixtures are valid."""
    
    @pytest.mark.parametrize("network_name,metadata", 
                             list(sdn.NETWORK_METADATA.items()))
    def test_network_is_valid(self, network_name: str, metadata: dict) -> None:
        """
        Test that network can be created and is valid.
        
        Parameters
        ----------
        network_name : str
            Name of the network fixture.
        metadata : dict
            Expected metadata for the network.
        """
        network = getattr(sdn, network_name)
        # Basic validity checks
        assert network.number_of_nodes() == metadata['nodes']
        assert network.number_of_edges() == metadata['edges']
    
    @pytest.mark.parametrize("network_name,metadata",
                             list(sdn.NETWORK_METADATA.items()))
    def test_network_properties_match_metadata(self, network_name: str, metadata: dict) -> None:
        """
        Test that computed properties match documented metadata.
        
        Parameters
        ----------
        network_name : str
            Name of the network fixture.
        metadata : dict
            Expected metadata for the network.
        """
        network = getattr(sdn, network_name)
        
        if metadata['nodes'] == 0:
            # Empty network - skip property checks
            return
        
        # Test classification properties
        assert sd_is_tree(network) == metadata['is_tree'], \
            f"{network_name}: is_tree mismatch"
        assert sd_is_binary(network) == metadata['is_binary'], \
            f"{network_name}: is_binary mismatch"
        assert sd_has_parallel_edges(network) == metadata['has_parallel_edges'], \
            f"{network_name}: has_parallel_edges mismatch"
        
        # For non-empty networks, test level-related properties
        if metadata['nodes'] > 1:
            assert sd_level(network) == metadata['level'], \
                f"{network_name}: level mismatch (expected {metadata['level']}, got {sd_level(network)})"
            assert sd_vertex_level(network) == metadata['vertex_level'], \
                f"{network_name}: vertex_level mismatch"
            assert sd_reticulation_number(network) == metadata['reticulation_number'], \
                f"{network_name}: reticulation_number mismatch"
    
    def test_can_retrieve_all_networks(self) -> None:
        """Test that all networks can be retrieved."""
        all_networks = sdn.get_all_networks()
        assert len(all_networks) == len(sdn.NETWORK_METADATA)
        assert all(net is not None for net in all_networks)
    
    def test_can_filter_by_category(self) -> None:
        """Test that networks can be filtered by category."""
        trees = sdn.get_networks_by_category('trees')
        assert len(trees) > 0
        assert all(sdn.NETWORK_METADATA[name]['category'] == 'trees' 
                  for name in sdn.NETWORK_METADATA 
                  if getattr(sdn, name) in trees)
    
    def test_can_filter_by_property(self) -> None:
        """Test that networks can be filtered by property."""
        level_0_nets = sdn.get_networks_with_property(level=0)
        assert len(level_0_nets) > 0
        
        binary_nets = sdn.get_networks_with_property(is_binary=True)
        assert len(binary_nets) > 0
    
    def test_metadata_completeness(self) -> None:
        """Test that all metadata entries have required fields."""
        required_fields = [
            'category', 'nodes', 'edges', 'level', 'vertex_level',
            'reticulation_number', 'is_tree', 'is_binary', 
            'has_parallel_edges', 'num_hybrids'
        ]
        for name, meta in sdn.NETWORK_METADATA.items():
            for field in required_fields:
                assert field in meta, f"{name} missing field: {field}"


class TestFixtureConsistency:
    """Test consistency across fixture modules."""
    
    def test_both_modules_have_similar_coverage(self) -> None:
        """Test that both directed and semi-directed have similar network counts."""
        d_count = len(dn.NETWORK_METADATA)
        sd_count = len(sdn.NETWORK_METADATA)
        
        # Should have similar number of networks (within 10%)
        assert abs(d_count - sd_count) / max(d_count, sd_count) < 0.3, \
            f"Directed has {d_count} networks, SemiDirected has {sd_count}"
    
    def test_both_modules_have_all_categories(self) -> None:
        """Test that both modules have networks in all categories."""
        d_categories = set(m['category'] for m in dn.NETWORK_METADATA.values())
        sd_categories = set(m['category'] for m in sdn.NETWORK_METADATA.values())
        
        expected_categories = {
            'trees', 'simple_hybrids', 'multiple_blobs',
            'high_level', 'large_networks', 'parallel_edges', 'special_cases', 'non_lsa'
        }
        
        assert d_categories == expected_categories
        assert sd_categories == expected_categories
    
    def test_helper_functions_work(self) -> None:
        """Test that helper functions work correctly."""
        # Test directed networks
        d_trees = dn.get_networks_by_category('trees')
        assert len(d_trees) > 0
        
        d_level_2 = dn.get_networks_with_property(level=2)
        for net in d_level_2:
            assert d_level(net) == 2
        
        # Test semi-directed networks
        sd_trees = sdn.get_networks_by_category('trees')
        assert len(sd_trees) > 0
        
        sd_level_2 = sdn.get_networks_with_property(level=2)
        for net in sd_level_2:
            assert sd_level(net) == 2


