"""
Tests for validation utilities.

Tests the validates decorator and no_validation_context.
"""

import pytest

from phylozoo.utils.validation import validates, no_validation_context


class TestValidatesDecorator:
    """Tests for the validates decorator."""

    def test_basic_validation(self):
        """Test that validation is called after mutation."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_count = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1
                if len(self._items) > 5:
                    raise ValueError("Too many items")

        obj = TestClass()
        
        # Add items - validation should be called
        obj.add_item(1)
        assert obj._validation_count == 1
        
        obj.add_item(2)
        assert obj._validation_count == 2
        
        # Add items up to limit
        obj.add_item(3)
        assert obj._validation_count == 3
        obj.add_item(4)
        assert obj._validation_count == 4
        obj.add_item(5)
        assert obj._validation_count == 5
        
        # Should raise error when limit exceeded
        with pytest.raises(ValueError, match="Too many items"):
            obj.add_item(6)

    def test_custom_validation_method_name(self):
        """Test that custom validation method names work."""
        class TestClass:
            def __init__(self):
                self._data = []
                self._validate_mode = True
                self._validation_depth = 0
                self._custom_validation_count = 0

            @validates('_custom_validate')
            def add_item(self, item):
                self._data.append(item)

            def _custom_validate(self):
                self._custom_validation_count += 1

        obj = TestClass()
        
        obj.add_item(1)
        assert obj._custom_validation_count == 1
        
        obj.add_item(2)
        assert obj._custom_validation_count == 2

    def test_validation_not_called_when_disabled(self):
        """Test that validation is not called when validation mode is disabled."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = False
                self._validation_depth = 0
                self._validation_count = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # Validation should not be called
        obj.add_item(1)
        obj.add_item(2)
        assert obj._validation_count == 0

    def test_validation_without_method(self):
        """Test that missing validation method doesn't cause error."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)
                # No _validate method defined

        obj = TestClass()
        
        # Should work fine - no validation method to call
        obj.add_item(1)
        obj.add_item(2)
        assert len(obj._items) == 2

    def test_validation_default_mode(self):
        """Test that validation mode defaults to True if not set."""
        class TestClass:
            def __init__(self):
                self._items = []
                # Don't set _validate_mode
                self._validation_count = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # Should validate by default
        obj.add_item(1)
        assert obj._validation_count == 1

    def test_validation_preserves_return_value(self):
        """Test that validation doesn't affect method return values."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)
                return len(self._items)

            def _validate(self):
                pass  # No-op validation

        obj = TestClass()
        
        result1 = obj.add_item(1)
        assert result1 == 1
        
        result2 = obj.add_item(2)
        assert result2 == 2

    def test_validation_with_exception(self):
        """Test that validation exceptions propagate correctly."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                if len(self._items) > 2:
                    raise ValueError("Limit exceeded")

        obj = TestClass()
        
        obj.add_item(1)  # OK
        obj.add_item(2)  # OK
        
        # Should raise error
        with pytest.raises(ValueError, match="Limit exceeded"):
            obj.add_item(3)
        
        # Item should still be added (validation happens after mutation)
        assert len(obj._items) == 3


class TestNoValidationContext:
    """Tests for the no_validation_context context manager."""

    def test_basic_no_validation(self):
        """Test that validation is disabled in context."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_count = 0

            def no_validation(self):
                return no_validation_context(self)

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # Normal operation - validation occurs
        obj.add_item(1)
        assert obj._validation_count == 1
        
        # In no_validation context - validation disabled
        with obj.no_validation():
            obj.add_item(2)
            obj.add_item(3)
            assert obj._validation_count == 1  # Not incremented
        
        # After context - validation should run once
        assert obj._validation_count == 2  # Called once on exit

    def test_nested_no_validation(self):
        """Test that nested no_validation contexts work correctly."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_count = 0

            def no_validation(self):
                return no_validation_context(self)

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # Nested contexts
        with obj.no_validation():
            obj.add_item(1)
            assert obj._validation_count == 0
            
            with obj.no_validation():
                obj.add_item(2)
                assert obj._validation_count == 0
            
            obj.add_item(3)
            assert obj._validation_count == 0
        
        # Validation should only run once when exiting outermost context
        assert obj._validation_count == 1

    def test_no_validation_initializes_state(self):
        """Test that no_validation initializes state if missing."""
        class TestClass:
            def __init__(self):
                # Don't initialize validation state
                self._items = []

            def no_validation(self):
                return no_validation_context(self)

            def _validate(self):
                pass

        obj = TestClass()
        
        # Should work even without initialization
        assert not hasattr(obj, '_validate_mode')
        assert not hasattr(obj, '_validation_depth')
        
        with obj.no_validation():
            pass
        
        # State should be initialized
        assert hasattr(obj, '_validate_mode')
        assert hasattr(obj, '_validation_depth')

    def test_no_validation_custom_method_name(self):
        """Test that custom validation method names work with no_validation."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._custom_validation_count = 0

            def no_validation(self):
                return no_validation_context(self, '_custom_validate')

            @validates('_custom_validate')
            def add_item(self, item):
                self._items.append(item)

            def _custom_validate(self):
                self._custom_validation_count += 1

        obj = TestClass()
        
        # Normal operation
        obj.add_item(1)
        assert obj._custom_validation_count == 1
        
        # In no_validation context
        with obj.no_validation():
            obj.add_item(2)
            assert obj._custom_validation_count == 1
        
        # Validation runs on exit
        assert obj._custom_validation_count == 2

    def test_no_validation_with_exception(self):
        """Test that exceptions in context still trigger validation on exit."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_count = 0

            def no_validation(self):
                return no_validation_context(self)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # Exception in context
        with pytest.raises(ValueError):
            with obj.no_validation():
                obj._items.append(1)
                raise ValueError("Test error")
        
        # Validation should still run on exit
        assert obj._validation_count == 1

    def test_no_validation_without_method(self):
        """Test that missing validation method doesn't cause error."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0

            def no_validation(self):
                return no_validation_context(self)
                # No _validate method defined

        obj = TestClass()
        
        # Should work fine
        with obj.no_validation():
            obj._items.append(1)

    def test_multiple_no_validation_contexts(self):
        """Test multiple separate no_validation contexts."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_count = 0

            def no_validation(self):
                return no_validation_context(self)

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # First context
        with obj.no_validation():
            obj.add_item(1)
            assert obj._validation_count == 0
        
        assert obj._validation_count == 1
        
        # Second context
        with obj.no_validation():
            obj.add_item(2)
            assert obj._validation_count == 1
        
        assert obj._validation_count == 2

    def test_no_validation_restores_mode(self):
        """Test that validation mode is restored after context."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0

            def no_validation(self):
                return no_validation_context(self)

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                pass

        obj = TestClass()
        
        assert obj._validate_mode is True
        
        with obj.no_validation():
            assert obj._validate_mode is False
        
        assert obj._validate_mode is True  # Restored

    def test_no_validation_with_initial_false_mode(self):
        """Test no_validation when validation mode starts as False."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = False
                self._validation_depth = 0
                self._validation_count = 0

            def no_validation(self):
                return no_validation_context(self)

            @validates()
            def add_item(self, item):
                self._items.append(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        # In context
        with obj.no_validation():
            obj.add_item(1)
            assert obj._validation_count == 0
        
        # Validation should be called on exit (no_validation_context calls it directly)
        # even though _validate_mode was False initially
        assert obj._validation_count == 1
        
        # Mode should be restored to False
        assert obj._validate_mode is False


class TestCombinedValidation:
    """Tests combining validation with other features."""

    def test_validation_with_multiple_mutations(self):
        """Test validation with multiple mutation methods."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_count = 0

            @validates()
            def add_item(self, item):
                self._items.append(item)

            @validates()
            def remove_item(self, item):
                if item in self._items:
                    self._items.remove(item)

            def _validate(self):
                self._validation_count += 1

        obj = TestClass()
        
        obj.add_item(1)
        assert obj._validation_count == 1
        
        obj.add_item(2)
        assert obj._validation_count == 2
        
        obj.remove_item(1)
        assert obj._validation_count == 3

    def test_validation_order(self):
        """Test that validation happens after mutation."""
        class TestClass:
            def __init__(self):
                self._items = []
                self._validate_mode = True
                self._validation_depth = 0
                self._validation_order = []

            @validates()
            def add_item(self, item):
                self._items.append(item)
                self._validation_order.append('mutation')

            def _validate(self):
                self._validation_order.append('validation')
                assert len(self._items) > 0  # Mutation already happened

        obj = TestClass()
        
        obj.add_item(1)
        assert obj._validation_order == ['mutation', 'validation']

