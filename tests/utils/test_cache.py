"""
Tests for caching utilities.

Tests the cached_property and clears_cache decorators.
"""

import pytest

from phylozoo.utils.cache import cached_property, clears_cache


class TestCachedProperty:
    """Tests for the cached_property decorator."""

    def test_basic_caching(self):
        """Test that cached properties are computed once and reused."""
        class TestClass:
            def __init__(self):
                self._data = [1, 2, 3]
                self._cache_version = 0
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return sum(self._data) ** 2

        obj = TestClass()
        
        # First access computes
        result1 = obj.expensive
        assert result1 == 36
        assert obj._compute_count == 1
        
        # Second access uses cache
        result2 = obj.expensive
        assert result2 == 36
        assert obj._compute_count == 1
        assert result1 == result2

    def test_cache_invalidation(self):
        """Test that cache is cleared when @clears_cache is called."""
        class TestClass:
            def __init__(self):
                self._data = [1, 2, 3]
                self._cache_version = 0
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return sum(self._data) ** 2

            @clears_cache
            def add_item(self, item):
                self._data.append(item)

        obj = TestClass()
        
        # First access
        result1 = obj.expensive
        assert result1 == 36
        assert obj._compute_count == 1
        
        # Modify data with @clears_cache
        obj.add_item(4)
        
        # Cache should be invalidated, recompute
        result2 = obj.expensive
        assert result2 == 100  # (1+2+3+4)^2
        assert obj._compute_count == 2

    def test_multiple_cached_properties(self):
        """Test that multiple cached properties work independently."""
        class TestClass:
            def __init__(self):
                self._data = [1, 2, 3]
                self._cache_version = 0
                self._sum_count = 0
                self._product_count = 0

            @cached_property
            def sum_squared(self):
                self._sum_count += 1
                return sum(self._data) ** 2

            @cached_property
            def product(self):
                self._product_count += 1
                result = 1
                for x in self._data:
                    result *= x
                return result

            @clears_cache
            def add_item(self, item):
                self._data.append(item)

        obj = TestClass()
        
        # Access both properties
        sum1 = obj.sum_squared
        prod1 = obj.product
        assert sum1 == 36
        assert prod1 == 6
        assert obj._sum_count == 1
        assert obj._product_count == 1
        
        # Access again - should use cache
        sum2 = obj.sum_squared
        prod2 = obj.product
        assert obj._sum_count == 1
        assert obj._product_count == 1
        
        # Modify data - both caches should be cleared
        obj.add_item(4)
        
        # Both should recompute
        sum3 = obj.sum_squared
        prod3 = obj.product
        assert sum3 == 100
        assert prod3 == 24
        assert obj._sum_count == 2
        assert obj._product_count == 2

    def test_mutable_types_return_copies(self):
        """Test that mutable types (set, list, dict) return copies."""
        class TestClass:
            def __init__(self):
                self._cache_version = 0

            @cached_property
            def cached_set(self):
                return {1, 2, 3}

            @cached_property
            def cached_list(self):
                return [1, 2, 3]

            @cached_property
            def cached_dict(self):
                return {'a': 1, 'b': 2}

        obj = TestClass()
        
        # Get cached values
        set1 = obj.cached_set
        list1 = obj.cached_list
        dict1 = obj.cached_dict
        
        # Modify the returned values
        set1.add(4)
        list1.append(4)
        dict1['c'] = 3
        
        # Get again - should be unchanged (we modified copies)
        set2 = obj.cached_set
        list2 = obj.cached_list
        dict2 = obj.cached_dict
        
        assert set2 == {1, 2, 3}
        assert list2 == [1, 2, 3]
        assert dict2 == {'a': 1, 'b': 2}
        
        # Objects should be different
        assert set1 is not set2
        assert list1 is not list2
        assert dict1 is not dict2

    def test_immutable_types_return_same(self):
        """Test that immutable types return the same object."""
        class TestClass:
            def __init__(self):
                self._cache_version = 0

            @cached_property
            def cached_int(self):
                return 42

            @cached_property
            def cached_str(self):
                return "hello"

            @cached_property
            def cached_tuple(self):
                return (1, 2, 3)

        obj = TestClass()
        
        # Get cached values
        int1 = obj.cached_int
        str1 = obj.cached_str
        tuple1 = obj.cached_tuple
        
        # Get again - should be same objects (immutable, so safe)
        int2 = obj.cached_int
        str2 = obj.cached_str
        tuple2 = obj.cached_tuple
        
        assert int1 is int2
        assert str1 is str2
        assert tuple1 is tuple2

    def test_cache_not_cleared_without_decorator(self):
        """Test that cache is not cleared when method lacks @clears_cache."""
        class TestClass:
            def __init__(self):
                self._data = [1, 2, 3]
                self._cache_version = 0
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return sum(self._data) ** 2

            def modify_data(self, item):
                """Method without @clears_cache."""
                self._data.append(item)

        obj = TestClass()
        
        # First access
        result1 = obj.expensive
        assert obj._compute_count == 1
        
        # Modify data without @clears_cache
        obj.modify_data(4)
        
        # Cache should still be valid
        result2 = obj.expensive
        assert result2 == 36  # Still old value (cache not cleared)
        assert obj._compute_count == 1  # Not recomputed

    def test_direct_attribute_modification(self):
        """Test that direct attribute modification doesn't clear cache."""
        class TestClass:
            def __init__(self):
                self._data = [1, 2, 3]
                self._cache_version = 0
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return sum(self._data) ** 2

        obj = TestClass()
        
        # First access
        result1 = obj.expensive
        assert obj._compute_count == 1
        
        # Direct modification
        obj._data.append(4)
        
        # Cache should still be valid (not cleared by direct modification)
        result2 = obj.expensive
        assert result2 == 36  # Still old value
        assert obj._compute_count == 1

    def test_multiple_clears_cache_calls(self):
        """Test that multiple @clears_cache calls increment version correctly."""
        class TestClass:
            def __init__(self):
                self._data = [1, 2, 3]
                self._cache_version = 0
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return sum(self._data) ** 2

            @clears_cache
            def add_item(self, item):
                self._data.append(item)

            @clears_cache
            def remove_item(self, item):
                if item in self._data:
                    self._data.remove(item)

        obj = TestClass()
        
        # First access
        result1 = obj.expensive
        assert obj._compute_count == 1
        assert obj._cache_version == 0
        
        # First mutation
        obj.add_item(4)
        assert obj._cache_version == 1
        
        # Access - should recompute
        result2 = obj.expensive
        assert obj._compute_count == 2
        
        # Second mutation
        obj.remove_item(4)
        assert obj._cache_version == 2
        
        # Access - should recompute again
        result3 = obj.expensive
        assert obj._compute_count == 3

    def test_cache_version_initialization(self):
        """Test that cache version is initialized if missing."""
        class TestClass:
            def __init__(self):
                # Don't initialize _cache_version
                self._data = [1, 2, 3]
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return sum(self._data) ** 2

        obj = TestClass()
        
        # Should work even without _cache_version in __init__
        assert not hasattr(obj, '_cache_version')
        result = obj.expensive
        assert result == 36
        assert hasattr(obj, '_cache_version')
        assert obj._cache_version == 0

    def test_clears_cache_initializes_version(self):
        """Test that @clears_cache initializes version if missing."""
        class TestClass:
            def __init__(self):
                # Don't initialize _cache_version
                self._data = [1, 2, 3]

            @clears_cache
            def add_item(self, item):
                self._data.append(item)

        obj = TestClass()
        
        # Should work even without _cache_version in __init__
        assert not hasattr(obj, '_cache_version')
        obj.add_item(4)
        assert hasattr(obj, '_cache_version')
        assert obj._cache_version == 1

    def test_empty_mutable_types(self):
        """Test that empty mutable types are handled correctly."""
        class TestClass:
            def __init__(self):
                self._cache_version = 0

            @cached_property
            def empty_set(self):
                return set()

            @cached_property
            def empty_list(self):
                return []

            @cached_property
            def empty_dict(self):
                return {}

        obj = TestClass()
        
        set1 = obj.empty_set
        list1 = obj.empty_list
        dict1 = obj.empty_dict
        
        # Modify
        set1.add(1)
        list1.append(1)
        dict1['a'] = 1
        
        # Get again - should still be empty (copies)
        set2 = obj.empty_set
        list2 = obj.empty_list
        dict2 = obj.empty_dict
        
        assert set2 == set()
        assert list2 == []
        assert dict2 == {}

    def test_nested_mutable_types(self):
        """Test that nested mutable types in tuples/lists are handled."""
        class TestClass:
            def __init__(self):
                self._cache_version = 0

            @cached_property
            def nested_list(self):
                return [[1, 2], [3, 4]]

            @cached_property
            def nested_dict(self):
                return {'a': {'x': 1}, 'b': {'y': 2}}

        obj = TestClass()
        
        list1 = obj.nested_list
        dict1 = obj.nested_dict
        
        # Modify nested structures
        list1[0].append(5)
        dict1['a']['z'] = 3
        
        # Get again - nested structures are still mutable (shallow copy)
        list2 = obj.nested_list
        dict2 = obj.nested_dict
        
        # Note: Only top-level is copied, nested structures are shared
        # This is expected behavior for shallow copying
        assert list2[0] == [1, 2, 5]  # Modified
        assert dict2['a'] == {'x': 1, 'z': 3}  # Modified

    def test_property_with_side_effects(self):
        """Test that cached properties with side effects only execute once."""
        class TestClass:
            def __init__(self):
                self._cache_version = 0
                self._side_effect_count = 0

            @cached_property
            def expensive_with_side_effect(self):
                self._side_effect_count += 1
                return self._side_effect_count

        obj = TestClass()
        
        # First access
        result1 = obj.expensive_with_side_effect
        assert result1 == 1
        assert obj._side_effect_count == 1
        
        # Second access - side effect should not occur again
        result2 = obj.expensive_with_side_effect
        assert result2 == 1  # Cached value
        assert obj._side_effect_count == 1  # Not incremented

    def test_concurrent_cache_access(self):
        """Test that cache works correctly with multiple rapid accesses."""
        class TestClass:
            def __init__(self):
                self._cache_version = 0
                self._compute_count = 0

            @cached_property
            def expensive(self):
                self._compute_count += 1
                return 42

        obj = TestClass()
        
        # Multiple rapid accesses
        results = [obj.expensive for _ in range(100)]
        
        # Should only compute once
        assert obj._compute_count == 1
        assert all(r == 42 for r in results)

