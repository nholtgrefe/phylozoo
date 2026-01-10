"""
Tests for cycle resolution module.
"""

import pytest

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.split.base import Split
from phylozoo.core.primitives.partition import Partition
from phylozoo.inference.squirrel.sqprofile import SqQuartetProfile
from phylozoo.inference.squirrel.sqprofileset import SqQuartetProfileSet
from phylozoo.inference.squirrel.cycle_resolution import _qprofiles_to_hybrid_ranking


class TestQProfilesToHybridRanking:
    """Tests for the _qprofiles_to_hybrid_ranking function."""

    def test_four_set_partition_with_cycle(self) -> None:
        """Test _qprofiles_to_hybrid_ranking with 4-set partition containing a cycle profile."""
        # Create a cycle profile (2 quartets) with reticulation_leaf
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = SqQuartetProfileSet(
            profiles=[(q1, 1.0), (q2, 1.0)],
            taxa={'A', 'B', 'C', 'D'}
        )
        # Get the profile and check it has 2 quartets
        profile = profileset.get_profile(frozenset({'A', 'B', 'C', 'D'}))
        assert profile is not None
        assert len(profile.quartets) == 2
        
        # Create profile with reticulation_leaf
        sq_profile = SqQuartetProfile(
            {q1: 1.0, q2: 1.0},
            reticulation_leaf='A'
        )
        profileset_with_ret = SqQuartetProfileSet(
            profiles=[sq_profile],
            taxa={'A', 'B', 'C', 'D'}
        )
        
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        order = _qprofiles_to_hybrid_ranking(profileset_with_ret, partition, weights=True)
        
        assert isinstance(order, list)
        assert len(order) == 4
        assert all(isinstance(s, frozenset) for s in order)
        # The set containing 'A' should be first (most likely reticulation)
        assert order[0] == frozenset({'A'})

    def test_four_set_partition_without_cycle(self) -> None:
        """Test _qprofile_to_hybrid_ranking with 4-set partition containing only split profiles."""
        # Create profiles with only 1 quartet each (splits, not cycles)
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profileset = SqQuartetProfileSet(
            profiles=[q1, q2],
            taxa={'A', 'B', 'C', 'D'}
        )
        
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)
        
        assert isinstance(order, list)
        assert len(order) == 4
        # All sets should have equal voting (0.0), order may be arbitrary
        assert all(isinstance(s, frozenset) for s in order)

    def test_larger_partition_with_cycles(self) -> None:
        """Test _qprofile_to_hybrid_ranking with partition larger than 4 sets."""
        # Create profiles for various 4-taxon combinations
        # Some with cycles (2 quartets), some with splits (1 quartet)
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))  # Split
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))  # Split
        q3 = Quartet(Split({'A', 'B'}, {'E', 'F'}))  # Split
        q4 = Quartet(Split({'A', 'E'}, {'B', 'F'}))  # Split
        
        # Create a cycle profile
        q5 = Quartet(Split({'C', 'D'}, {'E', 'F'}))
        q6 = Quartet(Split({'C', 'E'}, {'D', 'F'}))
        cycle_profile = SqQuartetProfile(
            {q5: 1.0, q6: 1.0},
            reticulation_leaf='C'
        )
        
        profileset = SqQuartetProfileSet(
            profiles=[q1, q2, q3, q4, cycle_profile],
            taxa={'A', 'B', 'C', 'D', 'E', 'F'}
        )
        
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)
        
        assert isinstance(order, list)
        assert len(order) == 6
        assert all(isinstance(s, frozenset) for s in order)
        # Sets involved in cycles should have higher voting scores
        # (exact order depends on cycle percentages)

    def test_weights_parameter(self) -> None:
        """Test _qprofile_to_hybrid_ranking with weights=True and weights=False."""
        # Create a cycle profile with weights
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        sq_profile = SqQuartetProfile(
            {q1: 2.0, q2: 2.0},
            reticulation_leaf='A'
        )
        profileset = SqQuartetProfileSet(
            profiles=[sq_profile],
            taxa={'A', 'B', 'C', 'D'}
        )
        
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        
        order_with_weights = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)
        order_without_weights = _qprofiles_to_hybrid_ranking(profileset, partition, weights=False)
        
        # Both should return the same order (A should be first)
        assert order_with_weights[0] == frozenset({'A'})
        assert order_without_weights[0] == frozenset({'A'})
        # But the voting scores would differ (not directly accessible)

    def test_empty_profileset(self) -> None:
        """Test _qprofile_to_hybrid_ranking with empty profile set."""
        profileset = SqQuartetProfileSet(profiles=[], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)
        
        assert isinstance(order, list)
        assert len(order) == 4
        # All sets should have equal voting (0.0)

    def test_partition_with_missing_profiles(self) -> None:
        """Test _qprofile_to_hybrid_ranking when some 4-taxon sets don't have profiles."""
        # Create profiles for only some 4-taxon combinations
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = SqQuartetProfileSet(
            profiles=[q1],
            taxa={'A', 'B', 'C', 'D', 'E', 'F'}
        )
        
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}, {'F'}])
        order = _qprofiles_to_hybrid_ranking(profileset, partition, weights=True)
        
        assert isinstance(order, list)
        assert len(order) == 6
        # Should handle missing profiles gracefully


class TestInsertCycle:
    """Tests for the _insert_cycle function."""

    def test_basic_insert_cycle(self) -> None:
        """Test _insert_cycle with a simple 3-leaf network."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from phylozoo.inference.squirrel.cycle_resolution import _insert_cycle
        
        # Create a simple network with cut-vertex
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        # Create circular set ordering
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        
        # Insert cycle with ranking
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)
        
        assert isinstance(new_net, SemiDirectedPhyNetwork)
        assert 3 not in new_net._graph.nodes  # Cut-vertex removed
        assert len(new_net._graph.nodes) > len(net._graph.nodes)  # New cycle nodes added
        assert new_net.number_of_edges() > net.number_of_edges()  # Cycle edges added

    def test_insert_cycle_with_four_leaves(self) -> None:
        """Test _insert_cycle with a 4-leaf network."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from phylozoo.inference.squirrel.cycle_resolution import _insert_cycle
        
        # Create a network with cut-vertex connecting 4 leaves
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 1), (5, 2), (5, 3), (5, 4)],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'})
            ]
        )
        
        # Create circular set ordering
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}, {'D'}])
        
        # Insert cycle with ranking
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'}), frozenset({'D'})]
        new_net = _insert_cycle(net, 5, ordering, reticulation_ranking=ranking)
        
        assert isinstance(new_net, SemiDirectedPhyNetwork)
        # Should have 4 cycle nodes (one per partition set)
        # Original nodes: 1, 2, 3, 4, 5 (cut-vertex)
        # New nodes should be: 1, 2, 3, 4 (leaves) + 4 cycle nodes
        original_leaf_nodes = {1, 2, 3, 4}
        new_nodes = set(new_net._graph.nodes)
        cycle_nodes = new_nodes - original_leaf_nodes
        assert len(cycle_nodes) == 4  # 4 cycle nodes
        assert 5 not in new_nodes or len(cycle_nodes) == 4  # Cut-vertex should be removed or replaced

    def test_insert_cycle_validation(self) -> None:
        """Test _insert_cycle validation errors."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from phylozoo.inference.squirrel.cycle_resolution import _insert_cycle
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'})]
        
        # Test non-cut-vertex error
        with pytest.raises(ValueError, match="is not a cut-vertex"):
            _insert_cycle(net, 1, ordering, reticulation_ranking=ranking)
        
        # Test mismatched ordering error
        wrong_ordering = CircularSetOrdering([{'A'}, {'B'}])
        with pytest.raises(ValueError, match="does not match partition"):
            _insert_cycle(net, 3, wrong_ordering, reticulation_ranking=ranking)
        
        # Test reticulation_set not in ordering error
        bad_ranking = [frozenset({'X'})]
        with pytest.raises(ValueError, match="in ranking is not part of the circular set ordering"):
            _insert_cycle(net, 3, ordering, reticulation_ranking=bad_ranking)

    def test_insert_cycle_preserves_taxa(self) -> None:
        """Test that _insert_cycle preserves all taxa."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from phylozoo.inference.squirrel.cycle_resolution import _insert_cycle
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)
        
        # All original taxa should still be present
        assert new_net.taxa == net.taxa
        assert 'A' in new_net.taxa
        assert 'B' in new_net.taxa
        assert 'C' in new_net.taxa

    def test_insert_cycle_returns_semidirected_when_valid(self) -> None:
        """Test _insert_cycle returns SemiDirectedPhyNetwork when valid hybrid found."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from phylozoo.inference.squirrel.cycle_resolution import _insert_cycle
        
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
        
        # Should return SemiDirectedPhyNetwork when valid hybrid found
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)
        
        assert isinstance(new_net, SemiDirectedPhyNetwork)
        assert 'A' in new_net.taxa

    def test_insert_cycle_without_ranking(self) -> None:
        """Test _insert_cycle without reticulation_ranking returns MixedPhyNetwork."""
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
        from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
        from phylozoo.inference.squirrel.cycle_resolution import _insert_cycle
        
        # Create a simple network
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
        
        # Should return MixedPhyNetwork when no ranking provided
        new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=None)
        assert isinstance(new_net, MixedPhyNetwork)
        assert not isinstance(new_net, SemiDirectedPhyNetwork)

