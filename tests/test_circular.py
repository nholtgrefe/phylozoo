"""
Tests for circular ordering classes.

This module contains comprehensive tests for CircularSetOrdering and CircularOrdering.
"""

import pytest
from typing import List, Set

from phylozoo.core.circular import CircularSetOrdering, CircularOrdering


class TestCircularSetOrderingCreation:
    """Test creation and initialization of CircularSetOrdering."""
    
    def test_basic_creation(self) -> None:
        """Test basic creation of CircularSetOrdering."""
        cso = CircularSetOrdering([{1, 2}, {3}, {4}])
        assert len(cso) == 3
        assert cso.size() == 4
        assert {1, 2} in cso.setorder
        assert {3} in cso.setorder
        assert {4} in cso.setorder
    
    def test_empty_ordering(self) -> None:
        """Test creation of empty CircularSetOrdering."""
        cso = CircularSetOrdering([])
        assert len(cso) == 0
        assert cso.size() == 0
        assert list(cso) == []
    
    def test_single_set(self) -> None:
        """Test creation with single set."""
        cso = CircularSetOrdering([{1, 2, 3}])
        assert len(cso) == 1
        assert cso.size() == 3
        assert {1, 2, 3} in cso.setorder
    
    def test_two_sets(self) -> None:
        """Test creation with two sets."""
        cso = CircularSetOrdering([{1}, {2}])
        assert len(cso) == 2
        assert cso.size() == 2
    
    def test_overlapping_sets_raises_error(self) -> None:
        """Test that overlapping sets raise ValueError."""
        with pytest.raises(ValueError, match="sets overlap"):
            CircularSetOrdering([{1, 2}, {2, 3}])
    
    def test_large_ordering(self) -> None:
        """Test creation with many sets."""
        sets = [{i} for i in range(10)]
        cso = CircularSetOrdering(sets)
        assert len(cso) == 10
        assert cso.size() == 10


class TestCircularSetOrderingCanonicalForm:
    """Test canonical form computation for CircularSetOrdering."""
    
    def test_canonical_form_same_regardless_of_input_order(self) -> None:
        """Test that canonical form is same for different input orders."""
        cso1 = CircularSetOrdering([{3, 4}, {1, 2}, {5}])
        cso2 = CircularSetOrdering([{5}, {1, 2}, {3, 4}])
        cso3 = CircularSetOrdering([{1, 2}, {5}, {3, 4}])
        
        assert cso1._setorder == cso2._setorder == cso3._setorder
    
    def test_canonical_form_cyclic_permutation(self) -> None:
        """Test that cyclic permutations produce same canonical form."""
        cso1 = CircularSetOrdering([{1}, {2}, {3}, {4}])
        cso2 = CircularSetOrdering([{2}, {3}, {4}, {1}])
        cso3 = CircularSetOrdering([{3}, {4}, {1}, {2}])
        cso4 = CircularSetOrdering([{4}, {1}, {2}, {3}])
        
        assert cso1._setorder == cso2._setorder == cso3._setorder == cso4._setorder
    
    def test_canonical_form_reversal(self) -> None:
        """Test that reversal produces same canonical form."""
        cso1 = CircularSetOrdering([{1}, {2}, {3}, {4}])
        cso2 = CircularSetOrdering([{4}, {3}, {2}, {1}])
        
        assert cso1._setorder == cso2._setorder
    
    def test_canonical_form_mixed_operations(self) -> None:
        """Test canonical form with various rotations and reversals."""
        base = [{1, 2}, {3}, {4, 5}, {6}]
        variations = [
            base,
            [{3}, {4, 5}, {6}, {1, 2}],  # Rotation
            [{6}, {4, 5}, {3}, {1, 2}],  # Reversal
            [{4, 5}, {6}, {1, 2}, {3}],  # Another rotation
        ]
        
        csos = [CircularSetOrdering(v) for v in variations]
        canonical_forms = [cso._setorder for cso in csos]
        
        # All should have the same canonical form
        assert all(cf == canonical_forms[0] for cf in canonical_forms)


class TestCircularSetOrderingEquality:
    """Test equality operations for CircularSetOrdering."""
    
    def test_equality_same_order(self) -> None:
        """Test equality with same order."""
        cso1 = CircularSetOrdering([{1, 2}, {3}, {4}])
        cso2 = CircularSetOrdering([{1, 2}, {3}, {4}])
        assert cso1 == cso2
    
    def test_equality_cyclic_permutation(self) -> None:
        """Test equality with cyclic permutation."""
        cso1 = CircularSetOrdering([{1}, {2}, {3}])
        cso2 = CircularSetOrdering([{2}, {3}, {1}])
        assert cso1 == cso2
    
    def test_equality_reversal(self) -> None:
        """Test equality with reversal."""
        cso1 = CircularSetOrdering([{1}, {2}, {3}])
        cso2 = CircularSetOrdering([{3}, {2}, {1}])
        assert cso1 == cso2
    
    def test_inequality_different_sets(self) -> None:
        """Test inequality with different sets."""
        cso1 = CircularSetOrdering([{1}, {2}])
        cso2 = CircularSetOrdering([{1}, {3}])
        assert cso1 != cso2
    
    def test_inequality_different_sizes(self) -> None:
        """Test inequality with different sizes."""
        cso1 = CircularSetOrdering([{1}, {2}])
        cso2 = CircularSetOrdering([{1}, {2}, {3}])
        assert cso1 != cso2
    
    def test_inequality_not_circular_set_ordering(self) -> None:
        """Test inequality with non-CircularSetOrdering."""
        cso = CircularSetOrdering([{1}, {2}])
        assert cso != "not a CircularSetOrdering"
        assert cso != 42
        assert cso != [{1}, {2}]


class TestCircularSetOrderingHashing:
    """Test hashing operations for CircularSetOrdering."""
    
    def test_hash_consistency(self) -> None:
        """Test that hash is consistent."""
        cso = CircularSetOrdering([{1, 2}, {3}, {4}])
        assert hash(cso) == hash(cso)
    
    def test_hash_same_for_equivalent_orderings(self) -> None:
        """Test that equivalent orderings have same hash."""
        cso1 = CircularSetOrdering([{1}, {2}, {3}])
        cso2 = CircularSetOrdering([{2}, {3}, {1}])
        cso3 = CircularSetOrdering([{3}, {2}, {1}])
        
        assert hash(cso1) == hash(cso2) == hash(cso3)
    
    def test_hash_different_for_different_orderings(self) -> None:
        """Test that different orderings have different hashes."""
        cso1 = CircularSetOrdering([{1}, {2}])
        cso2 = CircularSetOrdering([{1}, {3}])
        
        # Note: hash collisions are possible but unlikely for these cases
        assert hash(cso1) != hash(cso2)
    
    def test_hashable_in_set(self) -> None:
        """Test that CircularSetOrdering can be used in sets."""
        cso1 = CircularSetOrdering([{1}, {2}, {3}])
        cso2 = CircularSetOrdering([{2}, {3}, {1}])
        cso3 = CircularSetOrdering([{1}, {2}])
        
        s = {cso1, cso2, cso3}
        assert len(s) == 2  # cso1 and cso2 are equal


class TestCircularSetOrderingImmutability:
    """Test immutability of CircularSetOrdering."""
    
    def test_cannot_modify_setorder(self) -> None:
        """Test that setorder cannot be modified."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        
        with pytest.raises(AttributeError, match="immutable"):
            cso._setorder = ({4},)
    
    def test_cannot_modify_parts(self) -> None:
        """Test that parts cannot be modified."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        
        with pytest.raises(AttributeError, match="immutable"):
            cso._parts = ({4},)
    
    def test_cannot_add_attributes(self) -> None:
        """Test that new attributes cannot be added."""
        cso = CircularSetOrdering([{1}, {2}])
        
        with pytest.raises(AttributeError, match="immutable"):
            cso.new_attribute = "value"


class TestCircularSetOrderingIteration:
    """Test iteration over CircularSetOrdering."""
    
    def test_iteration_deterministic(self) -> None:
        """Test that iteration order is deterministic."""
        cso1 = CircularSetOrdering([{3, 4}, {1, 2}, {5}])
        cso2 = CircularSetOrdering([{5}, {1, 2}, {3, 4}])
        
        assert list(cso1) == list(cso2)
    
    def test_iteration_yields_frozensets(self) -> None:
        """Test that iteration yields frozensets."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        items = list(cso)
        
        assert all(isinstance(item, frozenset) for item in items)
        assert frozenset({1, 2}) in items
        assert frozenset({3}) in items
    
    def test_iteration_order_matches_setorder(self) -> None:
        """Test that iteration order matches setorder."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        assert list(cso) == list(cso.setorder)


class TestCircularSetOrderingAreNeighbors:
    """Test are_neighbors method."""
    
    def test_are_neighbors_adjacent(self) -> None:
        """Test that adjacent sets are neighbors."""
        cso = CircularSetOrdering([{1, 2}, {3}, {4}])
        assert cso.are_neighbors({1, 2}, {3}) is True
        assert cso.are_neighbors({3}, {4}) is True
    
    def test_are_neighbors_wrap_around(self) -> None:
        """Test that wrap-around sets are neighbors."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        assert cso.are_neighbors({1}, {3}) is True
        assert cso.are_neighbors({3}, {1}) is True
    
    def test_are_neighbors_not_adjacent(self) -> None:
        """Test that non-adjacent sets are not neighbors."""
        cso = CircularSetOrdering([{1}, {2}, {3}, {4}])
        assert cso.are_neighbors({1}, {3}) is False
        assert cso.are_neighbors({2}, {4}) is False
    
    def test_are_neighbors_same_set_raises_error(self) -> None:
        """Test that checking same set raises ValueError."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        
        with pytest.raises(ValueError, match="Cannot check if a set is neighbour of itself"):
            cso.are_neighbors({1}, {1})
    
    def test_are_neighbors_missing_set_raises_error(self) -> None:
        """Test that missing set raises ValueError."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        
        with pytest.raises(ValueError, match="not found in setorder"):
            cso.are_neighbors({1}, {4})
        
        with pytest.raises(ValueError, match="not found in setorder"):
            cso.are_neighbors({4}, {1})


class TestCircularSetOrderingAreSingletons:
    """Test are_singletons method."""
    
    def test_are_singletons_true(self) -> None:
        """Test are_singletons returns True for singleton sets."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        assert cso.are_singletons() is True
    
    def test_are_singletons_false(self) -> None:
        """Test are_singletons returns False for non-singleton sets."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        assert cso.are_singletons() is False
    
    def test_are_singletons_mixed(self) -> None:
        """Test are_singletons with mixed set sizes."""
        cso = CircularSetOrdering([{1}, {2, 3}, {4}])
        assert cso.are_singletons() is False


class TestCircularSetOrderingToCircularOrdering:
    """Test to_circular_ordering method."""
    
    def test_to_circular_ordering_success(self) -> None:
        """Test successful conversion to CircularOrdering."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        co = cso.to_circular_ordering()
        
        assert isinstance(co, CircularOrdering)
        assert len(co) == 3
        assert 1 in co
        assert 2 in co
        assert 3 in co
    
    def test_to_circular_ordering_raises_error(self) -> None:
        """Test that conversion raises error for non-singletons."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        
        with pytest.raises(ValueError, match="Not all sets are singletons"):
            cso.to_circular_ordering()


class TestCircularSetOrderingSuborderings:
    """Test suborderings method."""
    
    def test_suborderings_size_2(self) -> None:
        """Test suborderings of size 2."""
        cso = CircularSetOrdering([{1}, {2}, {3}, {4}])
        subs = list(cso.suborderings(size=2))
        
        assert len(subs) == 6  # C(4,2) = 6 combinations (not consecutive)
        assert all(isinstance(sub, CircularSetOrdering) for sub in subs)
        assert all(len(sub) == 2 for sub in subs)
    
    def test_suborderings_size_3(self) -> None:
        """Test suborderings of size 3."""
        cso = CircularSetOrdering([{1}, {2}, {3}, {4}])
        subs = list(cso.suborderings(size=3))
        
        assert len(subs) == 4  # C(4,3) = 4
        assert all(len(sub) == 3 for sub in subs)
    
    def test_suborderings_size_equals_length(self) -> None:
        """Test suborderings when size equals length."""
        cso = CircularSetOrdering([{1}, {2}, {3}])
        subs = list(cso.suborderings(size=3))
        
        assert len(subs) == 1
        assert subs[0] == cso


class TestCircularSetOrderingRepresentativeOrderings:
    """Test representative_orderings method."""
    
    def test_representative_orderings_singletons(self) -> None:
        """Test representative orderings with singleton sets."""
        cso = CircularSetOrdering([{1}, {2}])
        reps = list(cso.representative_orderings())
        
        assert len(reps) == 1
        assert isinstance(reps[0], CircularOrdering)
        assert list(reps[0]) == [1, 2]
    
    def test_representative_orderings_multiple_elements(self) -> None:
        """Test representative orderings with multiple elements per set."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        reps = list(cso.representative_orderings())
        
        assert len(reps) == 2  # 2 choices from first set * 1 from second
        assert all(isinstance(rep, CircularOrdering) for rep in reps)
        assert all(len(rep) == 2 for rep in reps)
    
    def test_representative_orderings_larger(self) -> None:
        """Test representative orderings with larger sets."""
        cso = CircularSetOrdering([{1, 2}, {3, 4}])
        reps = list(cso.representative_orderings())
        
        assert len(reps) == 4  # 2 * 2 = 4


class TestCircularOrderingCreation:
    """Test creation and initialization of CircularOrdering."""
    
    def test_basic_creation(self) -> None:
        """Test basic creation of CircularOrdering."""
        co = CircularOrdering([1, 2, 3, 4])
        assert len(co) == 4
        assert co.size() == 4
        assert 1 in co
        assert 2 in co
        assert 3 in co
        assert 4 in co
    
    def test_empty_ordering(self) -> None:
        """Test creation of empty CircularOrdering."""
        co = CircularOrdering([])
        assert len(co) == 0
        assert co.size() == 0
        assert list(co) == []
    
    def test_single_element(self) -> None:
        """Test creation with single element."""
        co = CircularOrdering([42])
        assert len(co) == 1
        assert co.size() == 1
        assert 42 in co
    
    def test_two_elements(self) -> None:
        """Test creation with two elements."""
        co = CircularOrdering([1, 2])
        assert len(co) == 2
        assert co.size() == 2
    
    def test_duplicate_elements_raises_error(self) -> None:
        """Test that duplicate elements raise ValueError."""
        with pytest.raises(ValueError, match="must be unique"):
            CircularOrdering([1, 2, 1])
    
    def test_different_numeric_types(self) -> None:
        """Test creation with different numeric types (all comparable)."""
        # Different numeric types that are comparable
        co = CircularOrdering([1, 2.5, 3])
        assert len(co) == 3
        assert 1 in co
        assert 2.5 in co
        assert 3 in co
    
    def test_mixed_types_raises_error(self) -> None:
        """Test that truly mixed types (str and int) raise TypeError in canonical form."""
        # The canonical form computation requires comparable types
        # This is expected behavior
        with pytest.raises(TypeError):
            CircularOrdering([1, "a", 2.5])


class TestCircularOrderingCanonicalForm:
    """Test canonical form computation for CircularOrdering."""
    
    def test_canonical_form_same_regardless_of_input_order(self) -> None:
        """Test that canonical form is same for different input orders."""
        co1 = CircularOrdering([3, 1, 2])
        co2 = CircularOrdering([1, 2, 3])
        co3 = CircularOrdering([2, 3, 1])
        
        assert co1._order == co2._order == co3._order
    
    def test_canonical_form_cyclic_permutation(self) -> None:
        """Test that cyclic permutations produce same canonical form."""
        co1 = CircularOrdering([1, 2, 3, 4])
        co2 = CircularOrdering([2, 3, 4, 1])
        co3 = CircularOrdering([3, 4, 1, 2])
        co4 = CircularOrdering([4, 1, 2, 3])
        
        assert co1._order == co2._order == co3._order == co4._order
    
    def test_canonical_form_reversal(self) -> None:
        """Test that reversal produces same canonical form."""
        co1 = CircularOrdering([1, 2, 3, 4])
        co2 = CircularOrdering([4, 3, 2, 1])
        
        assert co1._order == co2._order


class TestCircularOrderingEquality:
    """Test equality operations for CircularOrdering."""
    
    def test_equality_same_order(self) -> None:
        """Test equality with same order."""
        co1 = CircularOrdering([1, 2, 3])
        co2 = CircularOrdering([1, 2, 3])
        assert co1 == co2
    
    def test_equality_cyclic_permutation(self) -> None:
        """Test equality with cyclic permutation."""
        co1 = CircularOrdering([1, 2, 3])
        co2 = CircularOrdering([2, 3, 1])
        assert co1 == co2
    
    def test_equality_reversal(self) -> None:
        """Test equality with reversal."""
        co1 = CircularOrdering([1, 2, 3])
        co2 = CircularOrdering([3, 2, 1])
        assert co1 == co2
    
    def test_inequality_different_elements(self) -> None:
        """Test inequality with different elements."""
        co1 = CircularOrdering([1, 2])
        co2 = CircularOrdering([1, 3])
        assert co1 != co2
    
    def test_inequality_different_sizes(self) -> None:
        """Test inequality with different sizes."""
        co1 = CircularOrdering([1, 2])
        co2 = CircularOrdering([1, 2, 3])
        assert co1 != co2
    
    def test_inequality_not_circular_ordering(self) -> None:
        """Test inequality with non-CircularOrdering."""
        co = CircularOrdering([1, 2])
        assert co != "not a CircularOrdering"
        assert co != 42
        assert co != [1, 2]


class TestCircularOrderingHashing:
    """Test hashing operations for CircularOrdering."""
    
    def test_hash_consistency(self) -> None:
        """Test that hash is consistent."""
        co = CircularOrdering([1, 2, 3])
        assert hash(co) == hash(co)
    
    def test_hash_same_for_equivalent_orderings(self) -> None:
        """Test that equivalent orderings have same hash."""
        co1 = CircularOrdering([1, 2, 3])
        co2 = CircularOrdering([2, 3, 1])
        co3 = CircularOrdering([3, 2, 1])
        
        assert hash(co1) == hash(co2) == hash(co3)
    
    def test_hashable_in_set(self) -> None:
        """Test that CircularOrdering can be used in sets."""
        co1 = CircularOrdering([1, 2, 3])
        co2 = CircularOrdering([2, 3, 1])
        co3 = CircularOrdering([1, 2])
        
        s = {co1, co2, co3}
        assert len(s) == 2  # co1 and co2 are equal


class TestCircularOrderingImmutability:
    """Test immutability of CircularOrdering."""
    
    def test_cannot_modify_order(self) -> None:
        """Test that order cannot be modified."""
        co = CircularOrdering([1, 2, 3])
        
        with pytest.raises(AttributeError, match="immutable"):
            co._order = (4, 5, 6)
    
    def test_cannot_add_attributes(self) -> None:
        """Test that new attributes cannot be added."""
        co = CircularOrdering([1, 2])
        
        with pytest.raises(AttributeError, match="immutable"):
            co.new_attribute = "value"


class TestCircularOrderingIteration:
    """Test iteration over CircularOrdering."""
    
    def test_iteration_deterministic(self) -> None:
        """Test that iteration order is deterministic."""
        co1 = CircularOrdering([3, 1, 2])
        co2 = CircularOrdering([1, 2, 3])
        
        assert list(co1) == list(co2)
    
    def test_iteration_yields_elements(self) -> None:
        """Test that iteration yields elements."""
        co = CircularOrdering([1, 2, 3])
        items = list(co)
        
        assert items == [1, 2, 3] or items == [1, 2, 3]  # Canonical form
        assert all(isinstance(item, int) for item in items)
    
    def test_iteration_order_matches_order(self) -> None:
        """Test that iteration order matches order property."""
        co = CircularOrdering([1, 2, 3])
        assert list(co) == list(co.order)


class TestCircularOrderingAreNeighbors:
    """Test are_neighbors method for CircularOrdering."""
    
    def test_are_neighbors_adjacent(self) -> None:
        """Test that adjacent elements are neighbors."""
        co = CircularOrdering([1, 2, 3, 4])
        assert co.are_neighbors(1, 2) is True
        assert co.are_neighbors(2, 3) is True
        assert co.are_neighbors(3, 4) is True
    
    def test_are_neighbors_wrap_around(self) -> None:
        """Test that wrap-around elements are neighbors."""
        co = CircularOrdering([1, 2, 3])
        assert co.are_neighbors(1, 3) is True
        assert co.are_neighbors(3, 1) is True
    
    def test_are_neighbors_not_adjacent(self) -> None:
        """Test that non-adjacent elements are not neighbors."""
        co = CircularOrdering([1, 2, 3, 4])
        assert co.are_neighbors(1, 3) is False
        assert co.are_neighbors(2, 4) is False
    
    def test_are_neighbors_same_element_raises_error(self) -> None:
        """Test that checking same element raises ValueError."""
        co = CircularOrdering([1, 2, 3])
        
        with pytest.raises(ValueError, match="Cannot check if a set is neighbour of itself"):
            co.are_neighbors(1, 1)
    
    def test_are_neighbors_missing_element_raises_error(self) -> None:
        """Test that missing element raises ValueError."""
        co = CircularOrdering([1, 2, 3])
        
        with pytest.raises(ValueError, match="not found in setorder"):
            co.are_neighbors(1, 4)


class TestCircularOrderingSuborderings:
    """Test suborderings method for CircularOrdering."""
    
    def test_suborderings_size_2(self) -> None:
        """Test suborderings of size 2."""
        co = CircularOrdering([1, 2, 3, 4])
        subs = list(co.suborderings(size=2))
        
        assert len(subs) == 6  # C(4,2) = 6
        assert all(isinstance(sub, CircularOrdering) for sub in subs)
        assert all(len(sub) == 2 for sub in subs)
    
    def test_suborderings_size_3(self) -> None:
        """Test suborderings of size 3."""
        co = CircularOrdering([1, 2, 3, 4])
        subs = list(co.suborderings(size=3))
        
        assert len(subs) == 4  # C(4,3) = 4
        assert all(len(sub) == 3 for sub in subs)


class TestCircularOrderingToCircularSetOrdering:
    """Test to_circular_setordering method."""
    
    def test_to_circular_setordering_default(self) -> None:
        """Test conversion with default mapping (singletons)."""
        co = CircularOrdering([1, 2, 3])
        cso = co.to_circular_setordering()
        
        assert isinstance(cso, CircularSetOrdering)
        assert len(cso) == 3
        assert {1} in cso.setorder
        assert {2} in cso.setorder
        assert {3} in cso.setorder
    
    def test_to_circular_setordering_custom_mapping(self) -> None:
        """Test conversion with custom mapping."""
        co = CircularOrdering([1, 2, 3])
        mapping = {1: {1, 10}, 2: {2}, 3: {3, 30}}
        cso = co.to_circular_setordering(mapping)
        
        assert isinstance(cso, CircularSetOrdering)
        assert len(cso) == 3
        assert {1, 10} in cso.setorder
        assert {2} in cso.setorder
        assert {3, 30} in cso.setorder


class TestCircularOrderingEdgeCases:
    """Test edge cases for CircularOrdering."""
    
    def test_single_element_neighbors(self) -> None:
        """Test neighbors with single element (should be neighbors of itself)."""
        co = CircularOrdering([42])
        # With one element, it's its own neighbor (wrap-around)
        # But the method should raise an error for same element
        with pytest.raises(ValueError):
            co.are_neighbors(42, 42)
    
    def test_two_elements_always_neighbors(self) -> None:
        """Test that two elements are always neighbors."""
        co = CircularOrdering([1, 2])
        assert co.are_neighbors(1, 2) is True
        assert co.are_neighbors(2, 1) is True


class TestCircularSetOrderingEdgeCases:
    """Test edge cases for CircularSetOrdering."""
    
    def test_single_set_neighbors(self) -> None:
        """Test neighbors with single set."""
        cso = CircularSetOrdering([{1, 2, 3}])
        # With one set, it's its own neighbor (wrap-around)
        # But the method should raise an error for same set
        with pytest.raises(ValueError):
            cso.are_neighbors({1, 2, 3}, {1, 2, 3})
    
    def test_two_sets_always_neighbors(self) -> None:
        """Test that two sets are always neighbors."""
        cso = CircularSetOrdering([{1}, {2}])
        assert cso.are_neighbors({1}, {2}) is True
        assert cso.are_neighbors({2}, {1}) is True


class TestCircularOrderingRepr:
    """Test string representation of CircularOrdering."""
    
    def test_repr_basic(self) -> None:
        """Test basic repr."""
        co = CircularOrdering([1, 2, 3])
        repr_str = repr(co)
        assert "CircularOrdering" in repr_str
        assert "1" in repr_str or "2" in repr_str or "3" in repr_str
    
    def test_repr_empty(self) -> None:
        """Test repr of empty ordering."""
        co = CircularOrdering([])
        assert repr(co) == "CircularOrdering([])"


class TestCircularSetOrderingRepr:
    """Test string representation of CircularSetOrdering."""
    
    def test_repr_basic(self) -> None:
        """Test basic repr."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        repr_str = repr(cso)
        assert "CircularSetOrdering" in repr_str
    
    def test_repr_empty(self) -> None:
        """Test repr of empty ordering."""
        cso = CircularSetOrdering([])
        assert repr(cso) == "CircularSetOrdering([])"


class TestCircularOrderingContains:
    """Test __contains__ for CircularOrdering."""
    
    def test_contains_true(self) -> None:
        """Test that elements are found."""
        co = CircularOrdering([1, 2, 3])
        assert 1 in co
        assert 2 in co
        assert 3 in co
    
    def test_contains_false(self) -> None:
        """Test that non-elements are not found."""
        co = CircularOrdering([1, 2, 3])
        assert 4 not in co
        assert "a" not in co


class TestCircularSetOrderingContains:
    """Test __contains__ for CircularSetOrdering."""
    
    def test_contains_true(self) -> None:
        """Test that sets are found."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        assert {1, 2} in cso
        assert {3} in cso
        assert frozenset({1, 2}) in cso
    
    def test_contains_false(self) -> None:
        """Test that non-sets are not found."""
        cso = CircularSetOrdering([{1, 2}, {3}])
        assert {4} not in cso
        assert {1, 3} not in cso

