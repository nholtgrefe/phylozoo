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
    with no_validation(only=["validate"]):
        obj.validate()
    assert obj.validate_calls == 1


def test_validate_not_suppressed_for_disallowed_pattern() -> None:
    """
    Non-allowed patterns do not suppress methods.
    """
    obj = Dummy()
    with no_validation(only=["not_validate"]):
        obj.validate()
    assert obj.validate_calls == 2


def test_only_allowed_patterns_are_applied() -> None:
    """
    Patterns passed to no_validation are ignored if not in allowed set.
    """
    obj = Dummy()
    with no_validation(only=["not_allowed", "validate*"]):
        obj.validate()
    assert obj.validate_calls == 1
    obj.validate()
    assert obj.validate_calls == 2


def test_nested_validation_calls_inner_suppressed_only() -> None:
    """
    Suppress only inner validation while outer still runs.
    """
    obj = MultiValidate()
    with no_validation(only=["_check"]):
        obj.validate()
    assert obj.validate_calls == 1
    assert obj.check_calls == 0


def test_nested_validation_calls_outer_suppressed() -> None:
    """
    Suppress outer validation (and thus inner call) via pattern match.
    """
    obj = MultiValidate()
    with no_validation(only=["validate"]):
        obj.validate()
    assert obj.validate_calls == 0
    assert obj.check_calls == 0


def test_nested_contexts_restore_state() -> None:
    """
    Nested contexts suppress validation and restore afterward.
    """
    obj = Dummy()
    with no_validation(only=["validate"]):
        obj.validate()
        with no_validation(only=["validate"]):
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

