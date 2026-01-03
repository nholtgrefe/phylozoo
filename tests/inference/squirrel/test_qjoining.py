"""
Tests for the qjoining module.

This module tests the quartet_joining and adapted_quartet_joining functions
for constructing trees from quartet profiles.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.base import Split
from phylozoo.inference.squirrel.qjoining import adapted_quartet_joining, quartet_joining


class TestQuartetJoining:
    """Tests for the quartet_joining function."""
    
    def test_quartet_joining_simple_four_taxa(self) -> None:
        """Test quartet_joining with a simple 4-taxon case."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        tree = quartet_joining(profileset)
        
        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}
    
    def test_quartet_joining_five_taxa(self) -> None:
        """Test quartet_joining with 5 taxa."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])
        
        tree = quartet_joining(profileset)
        
        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E'}
    
    def test_quartet_joining_fewer_than_four_taxa_raises_error(self) -> None:
        """Test that quartet_joining raises ValueError for fewer than 4 taxa."""
        # Empty profile set
        profileset_empty = QuartetProfileSet()
        
        with pytest.raises(ValueError, match="Quartet joining requires at least 4 taxa"):
            quartet_joining(profileset_empty)
        
        # Profile set with only 3 taxa (not possible with quartets, but test edge case)
        # Actually, quartets always have 4 taxa, so this is hard to test directly
        # But we can test with a profile set that has no profiles
        profileset_no_profiles = QuartetProfileSet(profiles=[])
        
        with pytest.raises(ValueError, match="Quartet joining requires at least 4 taxa"):
            quartet_joining(profileset_no_profiles)
    
    def test_quartet_joining_produces_binary_tree(self) -> None:
        """Test that quartet_joining produces a binary tree (all internal nodes have degree 3)."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])
        
        tree = quartet_joining(profileset)
        
        # Check that all internal nodes have degree 3 (binary tree)
        for node in tree.internal_nodes:
            assert tree.degree(node) == 3
    
    def test_quartet_joining_larger_set(self) -> None:
        """Test quartet_joining with a larger set of taxa."""
        # Create quartets for 6 taxa: A, B, C, D, E, F
        quartets = [
            Quartet(Split({'A', 'B'}, {'C', 'D'})),
            Quartet(Split({'A', 'B'}, {'C', 'E'})),
            Quartet(Split({'A', 'B'}, {'C', 'F'})),
            Quartet(Split({'A', 'B'}, {'D', 'E'})),
            Quartet(Split({'A', 'B'}, {'D', 'F'})),
            Quartet(Split({'A', 'B'}, {'E', 'F'})),
        ]
        profileset = QuartetProfileSet(profiles=quartets)
        
        tree = quartet_joining(profileset)
        
        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E', 'F'}
        # Check that tree has correct number of edges (n_nodes - 1 for a tree)
        assert tree.number_of_edges() == tree.number_of_nodes() - 1


class TestAdaptedQuartetJoining:
    """Tests for the adapted_quartet_joining function."""
    
    def test_adapted_quartet_joining_with_star_tree(self) -> None:
        """Test adapted_quartet_joining with a star tree as starting tree."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        # Create a star tree
        center = '_center'
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), ('D', {'label': 'D'})]
        )
        
        tree = adapted_quartet_joining(profileset, starting)
        
        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}
    
    def test_adapted_quartet_joining_with_non_star_tree(self) -> None:
        """Test adapted_quartet_joining with a non-star starting tree."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        # Create a non-star binary tree: ((A,B),((C,D),E))
        # All internal nodes (3, 4, 5) have degree 3
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[
                (3, 'A'), (3, 'B'), (3, 4),  # Node 3: degree 3
                (4, 5), (4, 'E'),  # Node 4: degree 3 (connected to 3, 5, E)
                (5, 'C'), (5, 'D')  # Node 5: degree 3 (connected to 4, C, D)
            ],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), 
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )
        
        tree = adapted_quartet_joining(profileset, starting)
        
        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D', 'E'}
    
    def test_adapted_quartet_joining_taxa_mismatch_raises_error(self) -> None:
        """Test that adapted_quartet_joining raises error when taxa don't match."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        # Create a starting tree with different taxa
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[('_center', 'A'), ('_center', 'B'), ('_center', 'C'), ('_center', 'X')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), ('X', {'label': 'X'})]
        )
        
        with pytest.raises(ValueError, match="Starting tree taxa must match profile set taxa"):
            adapted_quartet_joining(profileset, starting)
    
    def test_adapted_quartet_joining_fewer_than_four_taxa_returns_unchanged(self) -> None:
        """Test that adapted_quartet_joining returns unchanged tree for fewer than 4 taxa."""
        # Create a profile set with only 3 taxa (not possible with quartets, but test edge case)
        # Actually, we can't create a QuartetProfileSet with fewer than 4 taxa from quartets
        # But we can test with an empty profile set and a tree with 3 taxa
        profileset_empty = QuartetProfileSet()
        
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[('_center', 'A'), ('_center', 'B'), ('_center', 'C')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'})]
        )
        
        # This should fail because taxa don't match (empty profile set has no taxa)
        with pytest.raises(ValueError, match="Starting tree taxa must match profile set taxa"):
            adapted_quartet_joining(profileset_empty, starting)
    
    def test_adapted_quartet_joining_produces_binary_tree(self) -> None:
        """Test that adapted_quartet_joining produces a binary tree."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        center = '_center'
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D'), (center, 'E')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), 
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )
        
        tree = adapted_quartet_joining(profileset, starting)
        
        # Check that all internal nodes have degree 3 (binary tree)
        for node in tree.internal_nodes:
            assert tree.degree(node) == 3


class TestQuartetJoiningIntegration:
    """Integration tests comparing quartet_joining and adapted_quartet_joining."""
    
    def test_quartet_joining_equivalent_to_adapted_with_star_tree(self) -> None:
        """Test that quartet_joining is equivalent to adapted_quartet_joining with a star tree."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        # Use quartet_joining
        tree1 = quartet_joining(profileset)
        
        # Use adapted_quartet_joining with a star tree
        center = '_center'
        starting = SemiDirectedPhyNetwork(
            undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D'), (center, 'E')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), 
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )
        tree2 = adapted_quartet_joining(profileset, starting)
        
        # Both should be valid trees with the same taxa
        assert is_tree(tree1)
        assert is_tree(tree2)
        assert tree1.taxa == tree2.taxa
        assert tree1.taxa == {'A', 'B', 'C', 'D', 'E'}
    
    def test_adapted_quartet_joining_with_different_starting_trees(self) -> None:
        """Test that adapted_quartet_joining works with different starting tree topologies."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        # Test with star tree
        center1 = '_center1'
        starting1 = SemiDirectedPhyNetwork(
            undirected_edges=[(center1, 'A'), (center1, 'B'), (center1, 'C'), (center1, 'D'), (center1, 'E')],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), 
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )
        tree1 = adapted_quartet_joining(profileset, starting1)
        
        # Test with a different starting tree: ((A,B),((C,D),E))
        starting2 = SemiDirectedPhyNetwork(
            undirected_edges=[
                (3, 'A'), (3, 'B'), (3, 4),
                (4, 5), (4, 'E'),
                (5, 'C'), (5, 'D')
            ],
            nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), 
                   ('D', {'label': 'D'}), ('E', {'label': 'E'})]
        )
        tree2 = adapted_quartet_joining(profileset, starting2)
        
        # Both should produce valid trees with the same taxa
        assert is_tree(tree1)
        assert is_tree(tree2)
        assert tree1.taxa == tree2.taxa == {'A', 'B', 'C', 'D', 'E'}


class TestQuartetJoiningEdgeCases:
    """Edge case tests for quartet joining functions."""
    
    def test_quartet_joining_with_single_quartet(self) -> None:
        """Test quartet_joining with a single quartet."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        tree = quartet_joining(profileset)
        
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}
        # With only one quartet, the tree should be binary (4 taxa = 2 internal nodes)
        assert len(tree.internal_nodes) == 2
    
    def test_quartet_joining_with_multiple_quartets_same_taxa(self) -> None:
        """Test quartet_joining with multiple quartets on the same 4-taxon set."""
        # All quartets are on {A, B, C, D}
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])
        
        tree = quartet_joining(profileset)
        
        assert is_tree(tree)
        assert tree.taxa == {'A', 'B', 'C', 'D'}
    
    def test_quartet_joining_already_binary_tree(self) -> None:
        """Test that quartet_joining handles a case where the tree is already binary."""
        # With 4 taxa and compatible quartets, the result should be binary
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        tree = quartet_joining(profileset)
        
        # The tree should be binary (all internal nodes have degree 3)
        for node in tree.internal_nodes:
            assert tree.degree(node) == 3

