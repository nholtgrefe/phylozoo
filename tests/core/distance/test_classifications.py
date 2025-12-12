"""
Tests for distance matrix classification functions.
"""

import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.distance.classifications import (
    is_metric,
    is_pseudo_metric,
    is_kalmanson,
    has_zero_diagonal,
    is_nonnegative,
)


class TestIsMetric:
    """Test is_metric function."""
    
    def test_metric_matrix(self) -> None:
        """Test that a valid metric matrix returns True."""
        # Euclidean distance matrix (satisfies all metric properties)
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is True
    
    def test_non_metric_violates_triangle_inequality(self) -> None:
        """Test that matrix violating triangle inequality returns False."""
        # d(A,C) = 5 > d(A,B) + d(B,C) = 1 + 1 = 2
        matrix = np.array([
            [0, 1, 5],
            [1, 0, 1],
            [5, 1, 0]
        ])
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is False
    
    def test_non_metric_negative_distance(self) -> None:
        """Test that matrix with negative distances returns False."""
        matrix = np.array([
            [0, -1],
            [-1, 0]
        ])
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is False
    
    def test_non_metric_non_zero_diagonal(self) -> None:
        """Test that matrix with non-zero diagonal returns False."""
        matrix = np.array([
            [1, 1],
            [1, 1]
        ])
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is False
    
    def test_single_label_metric(self) -> None:
        """Test single label matrix (always metric)."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is True


class TestIsPseudoMetric:
    """Test is_pseudo_metric function."""
    
    def test_pseudo_metric_matrix(self) -> None:
        """Test that a valid pseudo-metric matrix returns True."""
        # Matrix satisfying non-negativity and triangle inequality
        matrix = np.array([
            [0.1, 1, 2],
            [1, 0.1, 1],
            [2, 1, 0.1]
        ])
        dm = DistanceMatrix(matrix)
        assert is_pseudo_metric(dm) is True
    
    def test_metric_is_pseudo_metric(self) -> None:
        """Test that a metric is also a pseudo-metric."""
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix)
        assert is_pseudo_metric(dm) is True
    
    def test_non_pseudo_metric_violates_triangle(self) -> None:
        """Test that matrix violating triangle inequality returns False."""
        matrix = np.array([
            [0, 1, 5],
            [1, 0, 1],
            [5, 1, 0]
        ])
        dm = DistanceMatrix(matrix)
        assert is_pseudo_metric(dm) is False
    
    def test_non_pseudo_metric_negative(self) -> None:
        """Test that matrix with negative distances returns False."""
        matrix = np.array([
            [0, -1],
            [-1, 0]
        ])
        dm = DistanceMatrix(matrix)
        assert is_pseudo_metric(dm) is False


class TestIsKalmanson:
    """Test is_kalmanson function."""
    
    def test_kalmanson_matrix(self) -> None:
        """Test that a Kalmanson matrix returns True."""
        # Circular distance matrix (Kalmanson)
        matrix = np.array([
            [0, 1, 2, 2, 1],
            [1, 0, 1, 2, 2],
            [2, 1, 0, 1, 2],
            [2, 2, 1, 0, 1],
            [1, 2, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        assert is_kalmanson(dm, ['A', 'B', 'C', 'D', 'E']) is True
    
    def test_kalmanson_wrong_order(self) -> None:
        """Test that wrong circular order returns False."""
        matrix = np.array([
            [0, 1, 2, 2, 1],
            [1, 0, 1, 2, 2],
            [2, 1, 0, 1, 2],
            [2, 2, 1, 0, 1],
            [1, 2, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        # Wrong order
        assert is_kalmanson(dm, ['A', 'C', 'B', 'D', 'E']) is False
    
    def test_kalmanson_incomplete_order_raises_error(self) -> None:
        """Test that incomplete circular order raises ValueError."""
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        with pytest.raises(ValueError, match="must contain all labels"):
            is_kalmanson(dm, ['A', 'B'])
    
    def test_kalmanson_non_pseudo_metric_raises_error(self) -> None:
        """Test that non-pseudo-metric matrix raises ValueError."""
        # Matrix violating triangle inequality
        matrix = np.array([
            [0, 1, 5],
            [1, 0, 1],
            [5, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        with pytest.raises(ValueError, match="must be pseudo-metric"):
            is_kalmanson(dm, ['A', 'B', 'C'])
    
    def test_small_kalmanson_matrix(self) -> None:
        """Test Kalmanson check on small matrix."""
        # For n < 4, Kalmanson conditions are trivially satisfied
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        assert is_kalmanson(dm, ['A', 'B', 'C']) is True


class TestClassificationEdgeCases:
    """Test edge cases for classification functions."""
    
    def test_zero_matrix_is_metric(self) -> None:
        """Test that zero matrix is a metric."""
        matrix = np.zeros((3, 3))
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is True
        assert is_pseudo_metric(dm) is True
    
    def test_identity_matrix_is_metric(self) -> None:
        """Test that identity-like matrix (0 on diagonal, 1 elsewhere) is metric."""
        matrix = np.ones((3, 3))
        np.fill_diagonal(matrix, 0)
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is True
    
    def test_large_metric_matrix(self) -> None:
        """Test classification on larger matrix."""
        n = 10
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                matrix[i, j] = abs(i - j)
        dm = DistanceMatrix(matrix)
        assert is_metric(dm) is True
        assert is_pseudo_metric(dm) is True

