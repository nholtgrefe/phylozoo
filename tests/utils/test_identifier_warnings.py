"""
Tests for identifier warnings utility module.

This module tests the warning functions for Python keyword identifiers and values.
"""

import warnings

import pytest

from phylozoo.utils.identifier_warnings import warn_on_keyword, warn_on_none_value


class TestWarnOnKeyword:
    """Test cases for warn_on_keyword function."""

    def test_keyword_identifier_warns(self) -> None:
        """Keywords used as identifiers should emit a warning."""
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword("for", "Identifier")

    def test_keyword_key_warns(self) -> None:
        """Keywords used as keys should emit a warning."""
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword("class", "Key")

    def test_keyword_name_warns(self) -> None:
        """Keywords used as names should emit a warning."""
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword("def", "Name")

    def test_non_keyword_no_warning(self) -> None:
        """Non-keywords should not emit a warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            warn_on_keyword("my_identifier", "Identifier")
            warn_on_keyword("my_key", "Key")
            warn_on_keyword("my_name", "Name")

    def test_none_as_keyword_warns(self) -> None:
        """None is a keyword and should warn."""
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(None, "Identifier")

    def test_true_as_keyword_warns(self) -> None:
        """True is a keyword and should warn."""
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(True, "Identifier")

    def test_false_as_keyword_warns(self) -> None:
        """False is a keyword and should warn."""
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(False, "Identifier")

    def test_numeric_value_no_warning(self) -> None:
        """Numeric values should not emit a warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            warn_on_keyword(123, "Identifier")
            warn_on_keyword(45.6, "Identifier")

    def test_custom_context_in_message(self) -> None:
        """Custom context should appear in warning message."""
        with pytest.warns(UserWarning, match="Custom context.*is a Python keyword"):
            warn_on_keyword("if", "Custom context")

    def test_literal_python_keywords_warn(self) -> None:
        """Literal Python keywords (None, True, False) should warn."""
        # Test with actual Python keyword literals, not strings
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(None, "Identifier")  # Literal None
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(True, "Identifier")  # Literal True
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(False, "Identifier")  # Literal False

    def test_string_none_vs_literal_none(self) -> None:
        """String 'None' should warn (it's a keyword), but this tests warn_on_keyword."""
        # String "None" is also a keyword, so it should warn
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword("None", "Identifier")  # String "None"
        # Literal None should also warn
        with pytest.warns(UserWarning, match="is a Python keyword"):
            warn_on_keyword(None, "Identifier")  # Literal None


class TestWarnOnNoneValue:
    """Test cases for warn_on_none_value function."""

    def test_none_value_warns(self) -> None:
        """None values should emit a warning."""
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            warn_on_none_value(None, "Attribute 'weight'")

    def test_non_none_value_no_warning(self) -> None:
        """Non-None values should not emit a warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            warn_on_none_value(5.0, "Attribute 'weight'")
            warn_on_none_value("value", "Parameter 'x'")
            warn_on_none_value(True, "Setting 'mode'")
            warn_on_none_value(False, "Option 'debug'")

    def test_zero_no_warning(self) -> None:
        """Zero should not emit a warning (not None)."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            warn_on_none_value(0, "Attribute 'count'")

    def test_empty_string_no_warning(self) -> None:
        """Empty string should not emit a warning (not None)."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            warn_on_none_value("", "Attribute 'name'")

    def test_custom_context_in_message(self) -> None:
        """Custom context should appear in warning message."""
        with pytest.warns(UserWarning, match="Custom context.*has value None"):
            warn_on_none_value(None, "Custom context")

    def test_attribute_context_format(self) -> None:
        """Attribute context format should work correctly."""
        with pytest.warns(UserWarning, match="Attribute 'weight'.*has value None"):
            warn_on_none_value(None, "Attribute 'weight'")

    def test_parameter_context_format(self) -> None:
        """Parameter context format should work correctly."""
        with pytest.warns(UserWarning, match="Parameter 'x'.*has value None"):
            warn_on_none_value(None, "Parameter 'x'")

    def test_multiple_none_values_warn(self) -> None:
        """Multiple None values should each emit a warning."""
        with pytest.warns(UserWarning, match="has value None.*Python keyword") as record:
            warn_on_none_value(None, "Attribute 'weight'")
            warn_on_none_value(None, "Attribute 'length'")
        # Should have two warnings
        assert len(record) == 2

    def test_literal_none_vs_string_none(self) -> None:
        """Literal None should warn, but string 'None' should not (it's not None)."""
        # Literal None should warn
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            warn_on_none_value(None, "Attribute 'weight'")  # Literal None
        
        # String "None" is NOT the same as literal None, so should not warn
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            warn_on_none_value("None", "Attribute 'weight'")  # String "None", not literal None

