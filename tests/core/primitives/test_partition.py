"""
Tests for the partition module.
"""

import pytest
from phylozoo.core.primitives.partition import Partition


class TestPartitionCreation:
    """Test cases for Partition creation and basic properties."""

    def test_basic_creation(self) -> None:
        """Test creating a valid partition."""
        partition = Partition([{1, 2}, {3, 4}])
        assert len(partition) == 2
        assert partition.size() == 4
        assert partition.elements == frozenset({1, 2, 3, 4})

    def test_single_part(self) -> None:
        """Test partition with a single part."""
        partition = Partition([{1, 2, 3}])
        assert len(partition) == 1
        assert partition.size() == 3
        assert partition.elements == frozenset({1, 2, 3})

    def test_singleton_parts(self) -> None:
        """Test partition with all singleton parts."""
        partition = Partition([{1}, {2}, {3}])
        assert len(partition) == 3
        assert partition.size() == 3
        assert partition.elements == frozenset({1, 2, 3})

    def test_empty_partition(self) -> None:
        """Test partition with no parts."""
        partition = Partition([])
        assert len(partition) == 0
        assert partition.size() == 0
        assert partition.elements == frozenset()

    def test_part_with_empty_set(self) -> None:
        """Test partition containing an empty set."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            partition = Partition([{1, 2}, set()])
        assert len(partition) == 2
        assert partition.size() == 2
        assert partition.elements == frozenset({1, 2})

    def test_large_partition(self) -> None:
        """Test partition with many parts."""
        parts = [{i} for i in range(10)]
        partition = Partition(parts)
        assert len(partition) == 10
        assert partition.size() == 10

    def test_string_elements(self) -> None:
        """Test partition with string elements."""
        partition = Partition([{"a", "b"}, {"c", "d"}])
        assert len(partition) == 2
        assert partition.size() == 4
        assert "a" in partition.elements

    def test_mixed_types(self) -> None:
        """Test partition with mixed element types."""
        partition = Partition([{1, "a"}, {2, "b"}])
        assert len(partition) == 2
        assert partition.size() == 4


class TestPartitionValidation:
    """Test cases for partition validation."""

    def test_overlapping_sets_raises_error(self) -> None:
        """Test that overlapping sets raise ValueError."""
        with pytest.raises(ValueError, match="Invalid partition: sets overlap"):
            Partition([{1, 2}, {2, 3}])

    def test_duplicate_elements_raises_error(self) -> None:
        """Test that duplicate elements in different parts raise error."""
        with pytest.raises(ValueError, match="Invalid partition: sets overlap"):
            Partition([{1, 2}, {3, 4}, {1, 5}])

    def test_complete_overlap_raises_error(self) -> None:
        """Test that completely overlapping sets raise error."""
        with pytest.raises(ValueError, match="Invalid partition: sets overlap"):
            Partition([{1, 2, 3}, {1, 2, 3}])

    def test_partial_overlap_raises_error(self) -> None:
        """Test that partially overlapping sets raise error."""
        with pytest.raises(ValueError, match="Invalid partition: sets overlap"):
            Partition([{1, 2, 3}, {3, 4, 5}])

    def test_valid_partition_passes(self) -> None:
        """Test that valid disjoint sets pass validation."""
        partition = Partition([{1, 2}, {3, 4}, {5, 6}])
        assert partition._is_valid() is True


class TestPartitionImmutability:
    """Test cases for partition immutability."""

    def test_cannot_modify_parts(self) -> None:
        """Test that parts cannot be modified."""
        partition = Partition([{1, 2}, {3, 4}])
        with pytest.raises(AttributeError, match="Cannot set attribute"):
            partition.parts = [{5, 6}]

    def test_cannot_modify_elements(self) -> None:
        """Test that elements cannot be modified."""
        partition = Partition([{1, 2}, {3, 4}])
        with pytest.raises(AttributeError, match="Cannot set attribute"):
            partition.elements = frozenset({5, 6})

    def test_cannot_modify_private_attributes(self) -> None:
        """Test that private attributes cannot be modified."""
        partition = Partition([{1, 2}, {3, 4}])
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            partition._parts = tuple([frozenset({5})])

    def test_cannot_add_new_attributes(self) -> None:
        """Test that new attributes cannot be added."""
        partition = Partition([{1, 2}, {3, 4}])
        with pytest.raises(AttributeError, match="Cannot set attribute"):
            partition.new_attribute = "test"

    def test_parts_property_is_readonly(self) -> None:
        """Test that parts property returns immutable tuple."""
        partition = Partition([{1, 2}, {3, 4}])
        parts = partition.parts
        assert isinstance(parts, tuple)
        # Verify it's the same object on repeated access
        assert partition.parts is parts

    def test_elements_property_is_readonly(self) -> None:
        """Test that elements property returns immutable frozenset."""
        partition = Partition([{1, 2}, {3, 4}])
        elements = partition.elements
        assert isinstance(elements, frozenset)
        # Verify it's the same object on repeated access
        assert partition.elements is elements


class TestPartitionIteration:
    """Test cases for partition iteration."""

    def test_iteration_order_deterministic(self) -> None:
        """Test that iteration order is deterministic regardless of input order."""
        p1 = Partition([{3, 4}, {1, 2}, {5}])
        p2 = Partition([{5}, {1, 2}, {3, 4}])
        p3 = Partition([{1, 2}, {5}, {3, 4}])

        # All should iterate in the same order
        iter1 = list(p1)
        iter2 = list(p2)
        iter3 = list(p3)

        assert iter1 == iter2 == iter3

    def test_iteration_sorted_by_size_then_elements(self) -> None:
        """Test that iteration is sorted by size, then by elements."""
        partition = Partition([{3, 4}, {1, 2}, {5}])
        parts = list(partition)

        # Should be sorted: {5} (size 1), {1, 2} (size 2, min=1), {3, 4} (size 2, min=3)
        assert parts[0] == frozenset({5})
        assert parts[1] == frozenset({1, 2})
        assert parts[2] == frozenset({3, 4})


class TestPartitionCanonicalOrder:
    """Test cases for canonical ordering of partitions."""

    def test_canonical_order_same_parts_different_input(self) -> None:
        """Test that partitions with same parts in different input order have same canonical order."""
        p1 = Partition([{3, 4}, {1, 2}, {5}])
        p2 = Partition([{5}, {1, 2}, {3, 4}])
        p3 = Partition([{1, 2}, {5}, {3, 4}])

        # All should have identical _parts (canonical form)
        assert p1._parts == p2._parts == p3._parts

    def test_canonical_order_parts_property(self) -> None:
        """Test that parts property returns canonical order."""
        p1 = Partition([{3, 4}, {1, 2}, {5}])
        p2 = Partition([{5}, {1, 2}, {3, 4}])

        # Parts should be identical (canonical)
        assert p1.parts == p2.parts

    def test_canonical_order_repr_consistency(self) -> None:
        """Test that __repr__ is consistent for equivalent partitions."""
        p1 = Partition([{3, 4}, {1, 2}, {5}])
        p2 = Partition([{5}, {1, 2}, {3, 4}])
        p3 = Partition([{1, 2}, {5}, {3, 4}])

        # All should have identical string representation
        assert repr(p1) == repr(p2) == repr(p3)

    def test_canonical_order_stored_sorted(self) -> None:
        """Test that parts are stored in sorted canonical order."""
        partition = Partition([{3, 4}, {1, 2}, {5}])
        parts = partition._parts

        # Should be sorted: {5} (size 1), {1, 2} (size 2, min=1), {3, 4} (size 2, min=3)
        assert parts[0] == frozenset({5})
        assert parts[1] == frozenset({1, 2})
        assert parts[2] == frozenset({3, 4})

    def test_canonical_order_iteration_matches_parts(self) -> None:
        """Test that iteration order matches the stored canonical parts order."""
        partition = Partition([{3, 4}, {1, 2}, {5}])
        iter_parts = list(partition)
        stored_parts = list(partition._parts)

        # Iteration should match stored order exactly
        assert iter_parts == stored_parts

    def test_canonical_order_size_based(self) -> None:
        """Test that canonical order prioritizes size first."""
        partition = Partition([{1, 2, 3}, {4}, {5, 6}])
        parts = list(partition._parts)

        # Should be sorted by size: {4} (size 1), {5, 6} (size 2), {1, 2, 3} (size 3)
        assert len(parts[0]) == 1
        assert len(parts[1]) == 2
        assert len(parts[2]) == 3

    def test_canonical_order_elements_secondary(self) -> None:
        """Test that for same-size parts, elements are sorted."""
        partition = Partition([{3, 4}, {1, 2}, {5, 6}])
        parts = list(partition._parts)

        # All have size 2, so sorted by elements: {1, 2}, {3, 4}, {5, 6}
        assert parts[0] == frozenset({1, 2})
        assert parts[1] == frozenset({3, 4})
        assert parts[2] == frozenset({5, 6})

    def test_canonical_order_mixed_types(self) -> None:
        """Test that canonical order works with mixed types."""
        p1 = Partition([{1, "a"}, {2, "b"}])
        p2 = Partition([{2, "b"}, {1, "a"}])

        # Should have same canonical order despite different input
        assert p1._parts == p2._parts
        assert repr(p1) == repr(p2)

    def test_canonical_order_empty_partition(self) -> None:
        """Test that empty partition has empty canonical order."""
        partition = Partition([])
        assert partition._parts == tuple()
        assert list(partition) == []

    def test_canonical_order_single_part(self) -> None:
        """Test canonical order for single part partition."""
        partition = Partition([{1, 2, 3}])
        assert len(partition._parts) == 1
        assert partition._parts[0] == frozenset({1, 2, 3})

    def test_canonical_order_hash_consistency(self) -> None:
        """Test that hash is consistent with canonical order."""
        p1 = Partition([{3, 4}, {1, 2}])
        p2 = Partition([{1, 2}, {3, 4}])

        # Same canonical order should give same hash
        assert p1._parts == p2._parts
        assert hash(p1) == hash(p2)

    def test_iteration_preserves_frozensets(self) -> None:
        """Test that iteration returns frozensets."""
        partition = Partition([{1, 2}, {3, 4}])
        for part in partition:
            assert isinstance(part, frozenset)

    def test_iteration_covers_all_parts(self) -> None:
        """Test that iteration covers all parts."""
        partition = Partition([{1}, {2}, {3}, {4}])
        parts = list(partition)
        assert len(parts) == 4
        assert all(isinstance(p, frozenset) for p in parts)


class TestPartitionMethods:
    """Test cases for Partition methods."""

    def test_get_part_existing_element(self) -> None:
        """Test getting part for existing element."""
        partition = Partition([{1, 2}, {3, 4}])
        part = partition.get_part(1)
        assert part == frozenset({1, 2})

    def test_get_part_returns_correct_part(self) -> None:
        """Test that get_part returns the correct part."""
        partition = Partition([{1, 2}, {3, 4}, {5, 6}])
        assert partition.get_part(3) == frozenset({3, 4})
        assert partition.get_part(5) == frozenset({5, 6})

    def test_get_part_nonexistent_element_raises_error(self) -> None:
        """Test that get_part raises error for nonexistent element."""
        partition = Partition([{1, 2}, {3, 4}])
        with pytest.raises(ValueError, match="Element 5 not found in partition"):
            partition.get_part(5)

    def test_get_part_empty_partition_raises_error(self) -> None:
        """Test that get_part raises error on empty partition."""
        partition = Partition([])
        with pytest.raises(ValueError, match="Element 1 not found in partition"):
            partition.get_part(1)

    def test_subpartitions_size(self) -> None:
        """Test generating subpartitions of specified size."""
        partition = Partition([{1}, {2}, {3}, {4}])
        subparts = list(partition.subpartitions(size=2))
        assert len(subparts) == 6  # C(4,2) = 6

    def test_subpartitions_all_parts(self) -> None:
        """Test that subpartitions of full size returns original partition."""
        partition = Partition([{1, 2}, {3, 4}])
        subparts = list(partition.subpartitions(size=2))
        assert len(subparts) == 1
        assert subparts[0] == partition

    def test_subpartitions_size_one(self) -> None:
        """Test subpartitions of size 1."""
        partition = Partition([{1}, {2}, {3}])
        subparts = list(partition.subpartitions(size=1))
        assert len(subparts) == 3
        assert all(len(p) == 1 for p in subparts)

    def test_subpartitions_size_zero(self) -> None:
        """Test subpartitions of size 0 returns partition with empty parts."""
        partition = Partition([{1}, {2}, {3}])
        subparts = list(partition.subpartitions(size=0))
        # Size 0 combinations returns one empty combination, which is an empty partition
        assert len(subparts) == 1
        assert subparts[0] == Partition([])

    def test_subpartitions_size_too_large(self) -> None:
        """Test subpartitions with size larger than partition."""
        partition = Partition([{1}, {2}])
        subparts = list(partition.subpartitions(size=5))
        assert len(subparts) == 0

    def test_representative_partitions(self) -> None:
        """Test generating representative partitions."""
        partition = Partition([{1, 2}, {3, 4}])
        reps = list(partition.representative_partitions())
        assert len(reps) == 4  # 2 choices * 2 choices = 4

    def test_representative_partitions_singletons(self) -> None:
        """Test that representative partitions contain only singletons."""
        partition = Partition([{1, 2}, {3, 4}])
        reps = list(partition.representative_partitions())
        for rep in reps:
            assert all(len(part) == 1 for part in rep)

    def test_representative_partitions_covers_all_combinations(self) -> None:
        """Test that representative partitions cover all combinations."""
        partition = Partition([{1, 2}, {3}])
        reps = list(partition.representative_partitions())
        assert len(reps) == 2  # 2 choices * 1 choice = 2
        # Should have [{1}, {3}] and [{2}, {3}]
        expected = [
            Partition([{1}, {3}]),
            Partition([{2}, {3}]),
        ]
        assert set(reps) == set(expected)

    def test_is_refinement_true(self) -> None:
        """Test that fine partition is refinement of coarse partition."""
        fine = Partition([{1}, {2}, {3}])
        coarse = Partition([{1, 2}, {3}])
        assert fine.is_refinement(coarse) is True

    def test_is_refinement_false(self) -> None:
        """Test that coarse partition is not refinement of fine partition."""
        fine = Partition([{1}, {2}, {3}])
        coarse = Partition([{1, 2}, {3}])
        assert coarse.is_refinement(fine) is False

    def test_is_refinement_equal(self) -> None:
        """Test that equal partitions are refinements of each other."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{3, 4}, {1, 2}])
        assert p1.is_refinement(p2) is True
        assert p2.is_refinement(p1) is True

    def test_is_refinement_different_elements_raises_error(self) -> None:
        """Test that is_refinement raises error for different elements."""
        p1 = Partition([{1, 2}, {3}])
        p2 = Partition([{4, 5}, {6}])
        with pytest.raises(ValueError, match="Other partition covers different elements"):
            p1.is_refinement(p2)

    def test_is_refinement_not_partition_raises_error(self) -> None:
        """Test that is_refinement raises error for non-Partition."""
        partition = Partition([{1, 2}, {3, 4}])
        with pytest.raises(ValueError, match="The argument must be an instance of Partition"):
            partition.is_refinement("not a partition")


class TestPartitionEqualityAndHashing:
    """Test cases for partition equality and hashing."""

    def test_equality_same_order(self) -> None:
        """Test that partitions with same parts in same order are equal."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{1, 2}, {3, 4}])
        assert p1 == p2

    def test_equality_different_order(self) -> None:
        """Test that partitions with same parts in different order are equal."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{3, 4}, {1, 2}])
        assert p1 == p2

    def test_equality_different_parts(self) -> None:
        """Test that partitions with different parts are not equal."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{1, 2}, {5, 6}])
        assert p1 != p2

    def test_equality_not_partition(self) -> None:
        """Test that partition is not equal to non-Partition."""
        partition = Partition([{1, 2}, {3, 4}])
        assert partition != "not a partition"
        assert partition != [{1, 2}, {3, 4}]

    def test_hash_same_order(self) -> None:
        """Test that partitions with same parts have same hash."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{1, 2}, {3, 4}])
        assert hash(p1) == hash(p2)

    def test_hash_different_order(self) -> None:
        """Test that partitions with same parts in different order have same hash."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{3, 4}, {1, 2}])
        assert hash(p1) == hash(p2)

    def test_hashable_in_set(self) -> None:
        """Test that partitions can be used in sets."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{3, 4}, {1, 2}])
        p3 = Partition([{1, 2}, {5, 6}])
        partition_set = {p1, p2, p3}
        assert len(partition_set) == 2  # p1 and p2 are equal

    def test_hashable_in_dict(self) -> None:
        """Test that partitions can be used as dictionary keys."""
        p1 = Partition([{1, 2}, {3, 4}])
        p2 = Partition([{3, 4}, {1, 2}])
        d = {p1: "value1", p2: "value2"}
        assert len(d) == 1  # p1 and p2 are equal, so only one key
        assert d[p1] == "value2"  # Last value wins


class TestPartitionContains:
    """Test cases for the __contains__ method."""

    def test_contains_existing_part(self) -> None:
        """Test that __contains__ returns True for existing part."""
        partition = Partition([{1, 2}, {3, 4}])
        assert {1, 2} in partition
        assert {3, 4} in partition

    def test_contains_nonexistent_part(self) -> None:
        """Test that __contains__ returns False for nonexistent part."""
        partition = Partition([{1, 2}, {3, 4}])
        assert {1, 3} not in partition
        assert {5, 6} not in partition

    def test_contains_with_frozenset(self) -> None:
        """Test that __contains__ works with frozensets."""
        partition = Partition([{1, 2}, {3, 4}])
        assert frozenset({1, 2}) in partition
        assert frozenset({3, 4}) in partition

    def test_contains_partial_match(self) -> None:
        """Test that __contains__ requires exact match."""
        partition = Partition([{1, 2}, {3, 4}])
        assert {1} not in partition  # Subset but not exact match
        assert {1, 2, 5} not in partition  # Superset but not exact match


class TestPartitionRepr:
    """Test cases for partition string representation."""

    def test_repr_contains_partition(self) -> None:
        """Test that repr contains 'Partition'."""
        partition = Partition([{1, 2}, {3, 4}])
        repr_str = repr(partition)
        assert "Partition" in repr_str

    def test_repr_shows_parts(self) -> None:
        """Test that repr shows the parts."""
        partition = Partition([{1, 2}, {3, 4}])
        repr_str = repr(partition)
        # Should show the parts (order may vary in string representation)
        assert "1" in repr_str or "2" in repr_str


class TestPartitionEdgeCases:
    """Test cases for edge cases and special scenarios."""

    def test_empty_sets_in_partition(self) -> None:
        """Test partition with multiple empty sets."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            partition = Partition([set(), set(), {1}])
        assert len(partition) == 3
        assert partition.size() == 1

    def test_single_element_partition(self) -> None:
        """Test partition with single element."""
        partition = Partition([{1}])
        assert len(partition) == 1
        assert partition.size() == 1
        assert partition.get_part(1) == frozenset({1})

    def test_large_numbers(self) -> None:
        """Test partition with large numbers."""
        partition = Partition([{1000, 2000}, {3000, 4000}])
        assert len(partition) == 2
        assert partition.size() == 4

    def test_negative_numbers(self) -> None:
        """Test partition with negative numbers."""
        partition = Partition([{-1, -2}, {-3, -4}])
        assert len(partition) == 2
        assert partition.size() == 4

    def test_zero_in_partition(self) -> None:
        """Test partition containing zero."""
        partition = Partition([{0, 1}, {2, 3}])
        assert 0 in partition.elements
        assert partition.get_part(0) == frozenset({0, 1})


class TestPartitionDocstringExamples:
    """Test examples from Partition docstring."""
    
    def test_docstring_example(self) -> None:
        """Test the example from Partition docstring."""
        partition = Partition([{1, 2}, {3, 4}, {5}])
        assert len(partition) == 3
        assert partition.size() == 5
        assert {1, 2} in partition
        assert partition.get_part(3) == frozenset({3, 4})
        # Check that parts is a read-only tuple in canonical form
        assert isinstance(partition.parts, tuple)
        assert frozenset({1, 2}) in partition.parts
        assert frozenset({3, 4}) in partition.parts
        assert frozenset({5}) in partition.parts

