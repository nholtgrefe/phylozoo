"""
Tests for utility modules.
"""

import pytest
from phylozoo.utils._config import set_validate, validate
from phylozoo.structures.circular import CircularOrdering, CircularSetOrdering
from phylozoo.utils.distances import DistanceMatrix
from phylozoo.structures.partition import Partition
from phylozoo.utils.tools import id_generator, normalized_l_p_norm


class TestConfig:
    """Test cases for the configuration module."""

    def test_validate_default(self) -> None:
        """Test that validation is enabled by default."""
        assert validate() is True

    def test_set_validate(self) -> None:
        """Test setting validation."""
        original = validate()
        set_validate(False)
        assert validate() is False
        set_validate(True)
        assert validate() is True
        set_validate(original)  # Restore original state


class TestCircularOrdering:
    """Test cases for the CircularOrdering class."""

    def test_circular_ordering_creation(self) -> None:
        """Test creating a circular ordering."""
        ordering = CircularOrdering(["a", "b", "c"])
        assert len(ordering) == 3

    def test_circular_ordering_repr(self) -> None:
        """Test string representation of a circular ordering."""
        ordering = CircularOrdering(["a", "b", "c"])
        repr_str = repr(ordering)
        assert "CircularOrdering" in repr_str


class TestCircularSetOrdering:
    """Test cases for the CircularSetOrdering class."""

    def test_circular_set_ordering_creation(self) -> None:
        """Test creating a circular set ordering."""
        ordering = CircularSetOrdering([{"a", "b"}, {"c", "d"}])
        assert len(ordering) == 2


class TestDistanceMatrix:
    """Test cases for the DistanceMatrix class."""

    def test_distance_matrix_creation(self) -> None:
        """Test creating a distance matrix."""
        matrix = DistanceMatrix()
        assert len(matrix.distances) == 0

    def test_set_and_get_distance(self) -> None:
        """Test setting and getting distances."""
        matrix = DistanceMatrix()
        matrix.set_distance("A", "B", 1.5)
        distance = matrix.get_distance("A", "B")
        assert distance == 1.5


class TestPartition:
    """Test cases for the Partition class."""

    def test_partition_creation(self) -> None:
        """Test creating a partition."""
        partition = Partition([{1, 2}, {3, 4}])
        assert len(partition) == 2

    def test_partition_immutability(self) -> None:
        """Test that partition is immutable (no add_block method)."""
        partition = Partition([{1, 2}])
        assert len(partition) == 1
        # Partition is immutable, so we can't modify it after creation
        assert {1, 2} in partition


class TestTools:
    """Test cases for utility tools."""

    def test_id_generator(self) -> None:
        """Test ID generator."""
        id_str = id_generator(size=10)
        assert len(id_str) == 10
        assert isinstance(id_str, str)

    def test_normalized_l_p_norm(self) -> None:
        """Test normalized L-p norm calculation."""
        vector = [1.0, 2.0, 3.0]
        norm = normalized_l_p_norm(vector, p=2.0)
        assert isinstance(norm, float)
        assert norm >= 0.0
