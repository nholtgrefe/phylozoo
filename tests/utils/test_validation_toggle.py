"""
Tests for validation suppression utilities.
"""

import pytest

from phylozoo.utils.validation import no_validation, validation_aware


@validation_aware(allowed=["validate"], default=["validate"])
class Dummy:
    """
    Simple class to track validation calls.

    Attributes
    ----------
    validate_calls : int
        Counter of how many times ``validate`` was invoked.
    """

    def __init__(self) -> None:
        self.validate_calls = 0
        self.validate()

    def validate(self) -> bool:
        """Increment call count and return True."""
        self.validate_calls += 1
        return True


@validation_aware(allowed=["validate", "_check"], default=[])
class MultiValidate:
    """
    Class where one validation method calls another.

    Attributes
    ----------
    validate_calls : int
        Count of validate invocations.
    check_calls : int
        Count of _check invocations.
    """

    def __init__(self) -> None:
        self.validate_calls = 0
        self.check_calls = 0

    def validate(self) -> bool:
        """Run outer validation and call inner check."""
        self.validate_calls += 1
        self._check()
        return True

    def _check(self) -> bool:
        """Inner validation routine."""
        self.check_calls += 1
        return True


def test_validate_runs_by_default() -> None:
    """
    Validation runs during initialization by default.
    """
    obj = Dummy()
    assert obj.validate_calls == 1
    obj.validate()
    assert obj.validate_calls == 2


def test_validate_suppressed_during_init() -> None:
    """
    Validation is suppressed for initialization inside ``no_validation`` context.
    """
    with no_validation():
        obj = Dummy()
    assert obj.validate_calls == 0


def test_validate_suppressed_explicit_only() -> None:
    """
    Validation is suppressed when pattern matches allowed set.
    """
    obj = Dummy()
    with no_validation(methods=["validate"]):
        obj.validate()
    assert obj.validate_calls == 1


def test_validate_not_suppressed_for_disallowed_pattern() -> None:
    """
    Non-allowed patterns do not suppress methods.
    """
    obj = Dummy()
    with no_validation(methods=["not_validate"]):
        obj.validate()
    assert obj.validate_calls == 2


def test_only_allowed_patterns_are_applied() -> None:
    """
    Patterns passed to no_validation are ignored if not in allowed set.
    """
    obj = Dummy()
    with no_validation(methods=["not_allowed", "validate*"]):
        obj.validate()
    assert obj.validate_calls == 1
    obj.validate()
    assert obj.validate_calls == 2


def test_nested_validation_calls_inner_suppressed_only() -> None:
    """
    Suppress only inner validation while outer still runs.
    """
    obj = MultiValidate()
    with no_validation(methods=["_check"]):
        obj.validate()
    assert obj.validate_calls == 1
    assert obj.check_calls == 0


def test_nested_validation_calls_outer_suppressed() -> None:
    """
    Suppress outer validation (and thus inner call) via pattern match.
    """
    obj = MultiValidate()
    with no_validation(methods=["validate"]):
        obj.validate()
    assert obj.validate_calls == 0
    assert obj.check_calls == 0


def test_nested_contexts_restore_state() -> None:
    """
    Nested contexts suppress validation and restore afterward.
    """
    obj = Dummy()
    with no_validation(methods=["validate"]):
        obj.validate()
        with no_validation(methods=["validate"]):
            obj.validate()
    obj.validate()
    assert obj.validate_calls == 2


def test_invalid_default_pattern_raises() -> None:
    """
    Decorator raises if default patterns are not subset of allowed.
    """
    with pytest.raises(ValueError):
        @validation_aware(allowed=["validate"], default=["other"])
        class _Bad:
            def validate(self) -> bool:
                return True


def test_cross_class_initialization_suppressed() -> None:
    """
    Validation is suppressed for classes initialized from other classes' methods.
    
    This test verifies that the validation suppression context applies to all
    @validation_aware classes within the context, regardless of which class
    instantiates them or where they're instantiated.
    """
    @validation_aware(allowed=["validate"], default=["validate"])
    class ClassA:
        """Class that creates instances of ClassB."""
        
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        
        def validate(self) -> None:
            """Validation method for ClassA."""
            self.validated = True
        
        def create_b(self) -> 'ClassB':
            """Create an instance of ClassB from within ClassA's method."""
            return ClassB()
    
    @validation_aware(allowed=["validate"], default=["validate"])
    class ClassB:
        """Class that is created from ClassA."""
        
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        
        def validate(self) -> None:
            """Validation method for ClassB."""
            self.validated = True
    
    # Test that validation runs normally outside context
    a1 = ClassA()
    assert a1.validated is True
    b1 = a1.create_b()
    assert b1.validated is True
    
    # Test that validation is suppressed for both classes within context
    with no_validation():
        a2 = ClassA()
        assert a2.validated is False  # ClassA's validation suppressed
        b2 = a2.create_b()  # Created from within ClassA's method
        assert b2.validated is False  # ClassB's validation also suppressed
    
    # Test that validation works again after context
    a3 = ClassA()
    assert a3.validated is True


def test_class_specific_suppression() -> None:
    """
    Validation can be suppressed for specific classes only.
    """
    @validation_aware(allowed=["validate"], default=["validate"])
    class ClassA:
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
    
    @validation_aware(allowed=["validate"], default=["validate"])
    class ClassB:
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
    
    # Suppress only ClassA
    with no_validation(classes=["ClassA"]):
        a = ClassA()
        b = ClassB()
        assert a.validated is False  # ClassA suppressed
        assert b.validated is True   # ClassB NOT suppressed
    
    # Suppress only ClassB
    with no_validation(classes=["ClassB"]):
        a = ClassA()
        b = ClassB()
        assert a.validated is True   # ClassA NOT suppressed
        assert b.validated is False  # ClassB suppressed


def test_class_wildcard_suppression() -> None:
    """
    Class wildcards work for pattern matching.
    """
    @validation_aware(allowed=["validate"], default=["validate"])
    class TestClassA:
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
    
    @validation_aware(allowed=["validate"], default=["validate"])
    class TestClassB:
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
    
    @validation_aware(allowed=["validate"], default=["validate"])
    class OtherClass:
        def __init__(self) -> None:
            self.validated = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
    
    # Suppress all classes matching "TestClass*"
    with no_validation(classes=["TestClass*"]):
        a = TestClassA()
        b = TestClassB()
        c = OtherClass()
        assert a.validated is False  # Matches pattern
        assert b.validated is False  # Matches pattern
        assert c.validated is True   # Doesn't match pattern


def test_class_and_method_suppression() -> None:
    """
    Both class and method patterns can be specified together.
    """
    @validation_aware(allowed=["validate", "_check"], default=["validate"])
    class ClassA:
        def __init__(self) -> None:
            self.validated = False
            self.checked = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
            self._check()
        def _check(self) -> None:
            self.checked = True
    
    @validation_aware(allowed=["validate", "_check"], default=["validate"])
    class ClassB:
        def __init__(self) -> None:
            self.validated = False
            self.checked = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
            self._check()
        def _check(self) -> None:
            self.checked = True
    
    # Suppress only ClassA's validate method
    with no_validation(classes=["ClassA"], methods=["validate"]):
        a = ClassA()
        b = ClassB()
        assert a.validated is False  # ClassA.validate suppressed
        assert a.checked is False   # ClassA._check not called (validate was suppressed)
        assert b.validated is True   # ClassB.validate NOT suppressed
        assert b.checked is True     # ClassB._check runs normally (called from validate)
    
    # Test that _check can still be called directly on ClassA
    a._check()
    assert a.checked is True  # _check works when called directly


def test_class_suppression_with_default_methods() -> None:
    """
    Class suppression uses class defaults when methods not specified.
    """
    @validation_aware(allowed=["validate", "_check"], default=["validate"])
    class ClassA:
        def __init__(self) -> None:
            self.validated = False
            self.checked = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
        def _check(self) -> None:
            self.checked = True
    
    @validation_aware(allowed=["validate", "_check"], default=[])
    class ClassB:
        def __init__(self) -> None:
            self.validated = False
            self.checked = False
            self.validate()
        def validate(self) -> None:
            self.validated = True
            self._check()
        def _check(self) -> None:
            self.checked = True
    
    # Suppress ClassA only (uses its default: validate)
    with no_validation(classes=["ClassA"]):
        a = ClassA()
        b = ClassB()
        assert a.validated is False  # ClassA.validate suppressed (default)
        assert a.checked is False    # ClassA._check not called (validate was suppressed)
        assert b.validated is True    # ClassB NOT suppressed
        assert b.checked is True      # ClassB runs normally
    
    # Test that _check can still be called directly on ClassA
    a._check()
    assert a.checked is True  # _check works when called directly

