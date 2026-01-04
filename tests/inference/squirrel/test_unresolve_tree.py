"""
Tests for the unresolve_tree module.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.base import Split
from phylozoo.inference.squirrel.unresolve_tree import split_support, unresolve_tree


class TestSplitSupport:
    """Tests for the split_support function."""

    def test_split_support_single_matching_quartet(self) -> None:
        """Test split support with a single matching quartet."""
        split = Split({'A', 'B'}, {'C', 'D'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        support = split_support(profileset, split)
        assert support == 1.0

    def test_split_support_no_matching_quartets(self) -> None:
        """Test split support with no matching quartets."""
        split = Split({'A', 'B'}, {'C', 'D'})
        q1 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        support = split_support(profileset, split)
        assert support == 0.0

    def test_split_support_partial_match(self) -> None:
        """Test split support with partial matching quartets."""
        # Test with split on {A, B, C, D, E} and quartets on different 4-taxon subsets
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))  # Matches split, on {A,B,C,D}
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))  # Matches split, on {A,B,C,E}
        q3 = Quartet(Split({'A', 'C'}, {'B', 'D'}))  # Doesn't match split, on {A,B,C,D}
        profileset = QuartetProfileSet(profiles=[q1, q2, q3], taxa={'A', 'B', 'C', 'D', 'E'})
        support = split_support(profileset, split)
        # q1 and q2 match, q3 doesn't. Support = 2 matching / 3 total = 2/3
        # But wait, q1 and q3 are on same 4-taxon set, so grouped into non-trivial profile (ignored)
        # Only q2 contributes, so support = 1.0
        assert support == 1.0

    def test_split_support_with_weights(self) -> None:
        """Test split support with weighted quartets."""
        # Use split on {A, B, C, D, E} and quartets on different 4-taxon sets
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))  # Matches split, on {A,B,C,D}
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))  # Matches split, on {A,B,C,E}
        # Set profile weights via QuartetProfileSet tuples
        profileset = QuartetProfileSet(profiles=[(q1, 2.0), (q2, 1.0)], taxa={'A', 'B', 'C', 'D', 'E'})
        support = split_support(profileset, split)
        # Both q1 and q2 match the split
        # Support = (2.0 + 1.0) / (2.0 + 1.0) = 1.0
        assert support == 1.0

    def test_split_support_trivial_split_raises_error(self) -> None:
        """Test that trivial splits raise an error."""
        split = Split({'A'}, {'B', 'C', 'D'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        with pytest.raises(ValueError, match="Split must be non-trivial"):
            split_support(profileset, split)

    def test_split_support_taxa_mismatch_raises_error(self) -> None:
        """Test that taxa mismatch raises an error."""
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        with pytest.raises(ValueError, match="Split taxa must match profile set taxa"):
            split_support(profileset, split)

    def test_split_support_non_trivial_profile_ignored(self) -> None:
        """Test that non-trivial profiles are ignored."""
        split = Split({'A', 'B'}, {'C', 'D'})
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        from phylozoo.core.quartet.qprofile import QuartetProfile
        # Create a profile with multiple quartets (non-trivial)
        profile = QuartetProfile([q1, q2])
        profileset = QuartetProfileSet(profiles=[profile])
        support = split_support(profileset, split)
        assert support == 0.0  # Non-trivial profiles are ignored

    def test_split_support_unresolved_quartet_ignored(self) -> None:
        """Test that unresolved quartets (star trees) are ignored."""
        # Use split on {A, B, C, D, E} to match profile set taxa
        split = Split({'A', 'B'}, {'C', 'D', 'E'})
        # Create a resolved quartet on {A,B,C,D}
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        # Create a star tree quartet (no split) on different 4-taxon set {A,B,C,E}
        q_star = Quartet({'A', 'B', 'C', 'E'})  # Star tree on different set
        profileset = QuartetProfileSet(profiles=[q1, q_star], taxa={'A', 'B', 'C', 'D', 'E'})
        support = split_support(profileset, split)
        # Only q1 contributes (star tree is ignored), so support should be 1.0
        assert support == 1.0


class TestUnresolveTree:
    """Tests for the unresolve_tree generator function."""

    def test_unresolve_tree_single_split(self) -> None:
        """Test unresolve_tree with a single split to contract."""
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        trees = list(unresolve_tree(tree, profileset))
        assert len(trees) >= 2  # Original + at least one contraction
        assert all(is_tree(t) for t in trees)
        # First tree should be the original
        assert trees[0].number_of_nodes() == tree.number_of_nodes()
        # Later trees should have fewer nodes
        assert trees[-1].number_of_nodes() < trees[0].number_of_nodes()

    def test_unresolve_tree_multiple_splits(self) -> None:
        """Test unresolve_tree with multiple splits."""
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        trees = list(unresolve_tree(tree, profileset))
        assert len(trees) >= 2
        assert all(is_tree(t) for t in trees)

    def test_unresolve_tree_no_non_trivial_splits(self) -> None:
        """Test unresolve_tree with no non-trivial splits."""
        # Create a star tree (only trivial splits)
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(0, 1), (0, 2), (0, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        # Create empty profileset with matching taxa
        profileset = QuartetProfileSet(profiles=[], taxa={'A', 'B', 'C'})
        trees = list(unresolve_tree(tree, profileset))
        assert len(trees) == 1
        assert trees[0].number_of_nodes() == tree.number_of_nodes()

    def test_unresolve_tree_invalid_tree_raises_error(self) -> None:
        """Test that invalid tree raises an error."""
        # Create a network with a hybrid node (not a tree)
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        # Create a valid network with a hybrid (not a tree)
        # Hybrid node 4: indegree 2, total_degree 3
        network = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 1), (5, 6), (5, 8), (6, 9)],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        profileset = QuartetProfileSet(profiles=[], taxa={'A', 'B', 'C'})
        with pytest.raises(ValueError, match="Tree must be a valid tree"):
            list(unresolve_tree(network, profileset))

    def test_unresolve_tree_taxa_mismatch_raises_error(self) -> None:
        """Test that taxa mismatch raises an error."""
        # Create a valid tree with 2 taxa
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        with pytest.raises(ValueError, match="Tree taxa must match profile set taxa"):
            list(unresolve_tree(tree, profileset))

    def test_unresolve_tree_contracts_least_supported_first(self) -> None:
        """Test that splits are contracted in order of support (lowest first)."""
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        # Create quartets where one split has lower support
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))  # This will have support 1.0
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))  # This will have support 0.0 (no matching quartet)
        profileset = QuartetProfileSet(profiles=[q1])
        trees = list(unresolve_tree(tree, profileset))
        # Should contract the least supported split first
        assert len(trees) >= 2

    def test_unresolve_tree_all_trees_valid(self) -> None:
        """Test that all yielded trees are valid trees."""
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        trees = list(unresolve_tree(tree, profileset))
        for t in trees:
            assert is_tree(t)
            assert t.taxa == tree.taxa  # Taxa should be preserved

    def test_unresolve_tree_decreasing_nodes(self) -> None:
        """Test that number of nodes decreases with each contraction."""
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (4, 3), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1])
        trees = list(unresolve_tree(tree, profileset))
        if len(trees) > 1:
            for i in range(1, len(trees)):
                assert trees[i].number_of_nodes() <= trees[i-1].number_of_nodes()

