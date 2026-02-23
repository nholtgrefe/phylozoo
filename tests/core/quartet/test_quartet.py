"""
Tests for Quartet class.
"""

import pytest

from phylozoo.core.quartet import Quartet
from phylozoo.core.split.base import Split


class TestQuartetInit:
    """Tests for Quartet initialization."""
    
    def test_init_from_split(self) -> None:
        """Test creating a quartet from a split."""
        split = Split({1, 2}, {3, 4})
        quartet = Quartet(split)
        
        assert quartet.taxa == frozenset({1, 2, 3, 4})
        assert quartet.split == split
        assert quartet.is_resolved()
        assert not quartet.is_star()
    
    @pytest.mark.parametrize("taxa_input", [{1, 2, 3, 4}, frozenset({1, 2, 3, 4})])
    def test_init_from_set_like_star(self, taxa_input: set[int] | frozenset[int]) -> None:
        """Test creating a star quartet from set or frozenset of taxa."""
        quartet = Quartet(taxa_input)
        assert quartet.taxa == frozenset({1, 2, 3, 4})
        assert quartet.split is None
        assert not quartet.is_resolved()
        assert quartet.is_star()
    
    def test_init_trivial_split_error(self) -> None:
        """Test that trivial splits raise ValueError."""
        trivial_split = Split({1}, {2, 3, 4})
        
        with pytest.raises(ValueError, match="Split must be a 2\\|2 split"):
            Quartet(trivial_split)
    
    def test_init_wrong_number_taxa_split(self) -> None:
        """Test that splits with wrong number of taxa raise ValueError."""
        wrong_split = Split({1, 2}, {3, 4, 5})
        
        with pytest.raises(ValueError, match="Split must have exactly 4 elements"):
            Quartet(wrong_split)
    
    @pytest.mark.parametrize("taxa", [{1, 2, 3}, {1, 2, 3, 4, 5}])
    def test_init_wrong_number_taxa_star(self, taxa: set[int]) -> None:
        """Test that star trees with wrong number of taxa raise ValueError."""
        with pytest.raises(ValueError, match="Taxa must have exactly 4 elements"):
            Quartet(taxa)


class TestQuartetProperties:
    """Tests for Quartet properties."""
    
    def test_taxa_property(self) -> None:
        """Test taxa property."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        assert quartet.taxa == frozenset({1, 2, 3, 4})
        assert isinstance(quartet.taxa, frozenset)
    
    def test_split_property(self) -> None:
        """Test split property."""
        split = Split({1, 2}, {3, 4})
        quartet = Quartet(split)
        assert quartet.split == split
        
        star_quartet = Quartet({1, 2, 3, 4})
        assert star_quartet.split is None
    
    def test_immutability(self) -> None:
        """Test that quartet is immutable."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            quartet._taxa = frozenset({5, 6, 7, 8})
        
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            quartet._split = None


class TestQuartetMethods:
    """Tests for Quartet methods."""
    
    @pytest.mark.parametrize(
        "quartet_input,expected_resolved",
        [
            (Split({1, 2}, {3, 4}), True),
            ({1, 2, 3, 4}, False),
        ],
    )
    def test_is_resolved(
        self, quartet_input: Split | set[int], expected_resolved: bool
    ) -> None:
        """Test is_resolved method."""
        quartet = Quartet(quartet_input)
        assert quartet.is_resolved() == expected_resolved

    @pytest.mark.parametrize(
        "quartet_input,expected_star",
        [
            (Split({1, 2}, {3, 4}), False),
            ({1, 2, 3, 4}, True),
        ],
    )
    def test_is_star(
        self, quartet_input: Split | set[int], expected_star: bool
    ) -> None:
        """Test is_star method."""
        quartet = Quartet(quartet_input)
        assert quartet.is_star() == expected_star
    
    @pytest.mark.parametrize(
        "quartet_input",
        [Split({1, 2}, {3, 4}), {1, 2, 3, 4}],
    )
    def test_copy(self, quartet_input: Split | set[int]) -> None:
        """Test copy method for resolved and star quartets."""
        original = Quartet(quartet_input)
        copied = original.copy()
        assert copied is not original
        assert copied.taxa == original.taxa
        assert copied.split == original.split
        assert copied == original
    
    def test_to_network_resolved(self) -> None:
        """Test to_network for resolved quartet."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        network = quartet.to_network()
        
        assert network.number_of_nodes() == 6  # 4 leaves + 2 internal nodes
        assert network.number_of_edges() == 5  # 4 leaf edges + 1 connecting internal nodes
        assert network.taxa == {'1', '2', '3', '4'}
        assert len(network.leaves) == 4
        assert len(network.tree_nodes) == 2
    
    def test_to_network_star(self) -> None:
        """Test to_network for star tree."""
        quartet = Quartet({1, 2, 3, 4})
        network = quartet.to_network()
        
        assert network.number_of_nodes() == 5  # 4 leaves + 1 internal node
        assert network.number_of_edges() == 4  # 4 leaf edges
        assert network.taxa == {'1', '2', '3', '4'}
        assert len(network.leaves) == 4
        assert len(network.tree_nodes) == 1


class TestQuartetEquality:
    """Tests for Quartet equality and hashing."""
    
    def test_eq_same_split(self) -> None:
        """Test equality with same split."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 2}, {3, 4}))
        
        assert q1 == q2
        assert hash(q1) == hash(q2)
    
    def test_eq_different_split_same_taxa(self) -> None:
        """Test equality with different split but same taxa."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        assert q1 != q2
    
    def test_eq_same_star(self) -> None:
        """Test equality with same star tree."""
        q1 = Quartet({1, 2, 3, 4})
        q2 = Quartet({1, 2, 3, 4})
        
        assert q1 == q2
        assert hash(q1) == hash(q2)
    
    def test_eq_different_taxa(self) -> None:
        """Test equality with different taxa."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        assert q1 != q2
    
    def test_eq_resolved_vs_star(self) -> None:
        """Test equality between resolved and star on same taxa."""
        resolved = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        
        assert resolved != star
    
    def test_eq_different_type(self) -> None:
        """Test equality with different type."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        
        assert quartet != "not a quartet"
        assert quartet != 42
        assert quartet != Split({1, 2}, {3, 4})


class TestQuartetRepr:
    """Tests for Quartet string representation."""

    @pytest.mark.parametrize(
        "quartet_input,expect_split_in_repr,check_taxa_in_repr",
        [
            (Split({1, 2}, {3, 4}), True, False),
            ({1, 2, 3, 4}, False, True),
        ],
    )
    def test_repr(
        self,
        quartet_input: Split | set[int],
        expect_split_in_repr: bool,
        check_taxa_in_repr: bool,
    ) -> None:
        """Test string representation of resolved and star quartets."""
        quartet = Quartet(quartet_input)
        repr_str = repr(quartet)
        assert repr_str.startswith("Quartet(")
        if expect_split_in_repr:
            assert "Split" in repr_str
        else:
            assert "Split" not in repr_str
        if check_taxa_in_repr:
            assert "1" in repr_str and "2" in repr_str


class TestQuartetHashability:
    """Tests for Quartet hashability."""
    
    def test_hashable(self) -> None:
        """Test that quartets are hashable."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet({1, 2, 3, 4})
        
        # Should not raise
        hash(q1)
        hash(q2)
    
    def test_hash_consistency(self) -> None:
        """Test that hash is consistent."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 2}, {3, 4}))
        
        assert hash(q1) == hash(q2)
    
    def test_can_use_in_set(self) -> None:
        """Test that quartets can be used in sets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet({1, 2, 3, 4})
        
        quartet_set = {q1, q2, q3}
        assert len(quartet_set) == 3
        assert q1 in quartet_set
        assert q2 in quartet_set
        assert q3 in quartet_set
    
    def test_can_use_in_dict(self) -> None:
        """Test that quartets can be used as dictionary keys."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        quartet_dict = {q1: "first", q2: "second"}
        assert quartet_dict[q1] == "first"
        assert quartet_dict[q2] == "second"


class TestQuartetCircularOrderings:
    """Tests for Quartet circular_orderings property."""
    
    def test_circular_orderings_resolved(self) -> None:
        """Test circular_orderings for resolved quartet."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        orderings = quartet.circular_orderings
        
        assert len(orderings) == 2
        # Should contain [1,2,3,4] and [1,2,4,3] (where 1,2 are neighbors and 3,4 are neighbors)
        ordering_lists = [list(co.order) for co in orderings]
        assert [1, 2, 3, 4] in ordering_lists
        assert [1, 2, 4, 3] in ordering_lists
    
    def test_circular_orderings_resolved_different_split(self) -> None:
        """Test circular_orderings for resolved quartet with different split."""
        quartet = Quartet(Split({1, 3}, {2, 4}))
        orderings = quartet.circular_orderings
        
        assert len(orderings) == 2
        # Should contain orderings where 1,3 are neighbors and 2,4 are neighbors
        ordering_lists = [list(co.order) for co in orderings]
        # Check that 1 and 3 are neighbors in both orderings
        for ordering in orderings:
            assert ordering.are_neighbors(1, 3)
            assert ordering.are_neighbors(2, 4)
    
    def test_circular_orderings_star(self) -> None:
        """Test circular_orderings for star tree."""
        quartet = Quartet({1, 2, 3, 4})
        orderings = quartet.circular_orderings
        
        assert len(orderings) == 3
        # Should contain all three distinct orderings: abcd, acbd, abdc
        ordering_lists = [list(co.order) for co in orderings]
        assert [1, 2, 3, 4] in ordering_lists
        assert [1, 3, 2, 4] in ordering_lists
        assert [1, 2, 4, 3] in ordering_lists
    
    def test_circular_orderings_cached(self) -> None:
        """Test that circular_orderings is cached."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        
        # First access
        orderings1 = quartet.circular_orderings
        # Second access should return the same object (cached)
        orderings2 = quartet.circular_orderings
        
        assert orderings1 is orderings2
    
    def test_circular_orderings_resolved_neighbors(self) -> None:
        """Test that resolved quartet orderings have correct neighbor relationships."""
        quartet = Quartet(Split({1, 2}, {3, 4}))
        orderings = quartet.circular_orderings
        
        for ordering in orderings:
            # In both orderings, 1 and 2 should be neighbors
            assert ordering.are_neighbors(1, 2)
            # And 3 and 4 should be neighbors
            assert ordering.are_neighbors(3, 4)
    
    def test_circular_orderings_different_taxa(self) -> None:
        """Test circular_orderings with different taxa."""
        quartet = Quartet(Split({5, 6}, {7, 8}))
        orderings = quartet.circular_orderings
        
        assert len(orderings) == 2
        for ordering in orderings:
            assert ordering.are_neighbors(5, 6)
            assert ordering.are_neighbors(7, 8)
    
    def test_circular_orderings_star_all_orderings(self) -> None:
        """Test that star tree returns all three distinct circular orderings."""
        quartet = Quartet({1, 2, 3, 4})
        orderings = quartet.circular_orderings
        
        # Verify all three are distinct
        ordering_tuples = [tuple(co.order) for co in orderings]
        assert len(set(ordering_tuples)) == 3
        
        # Verify they are all valid circular orderings
        for ordering in orderings:
            assert len(ordering) == 4
            assert set(ordering) == {1, 2, 3, 4}

