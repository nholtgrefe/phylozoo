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
from phylozoo.core.primitives.circular_ordering import CircularOrdering
from phylozoo.utils.exceptions import PhyloZooValueError


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
        co = CircularOrdering(['A', 'B', 'C', 'D', 'E'])
        assert is_kalmanson(dm, co) is True
    
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
        co = CircularOrdering(['A', 'C', 'B', 'D', 'E'])
        assert is_kalmanson(dm, co) is False
    
    def test_kalmanson_incomplete_order_raises_error(self) -> None:
        """Test that incomplete circular order raises ValueError."""
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        co = CircularOrdering(['A', 'B'])
        with pytest.raises(PhyloZooValueError, match="must contain all labels"):
            is_kalmanson(dm, co)
    
    def test_kalmanson_non_pseudo_metric_raises_error(self) -> None:
        """Test that non-pseudo-metric matrix raises PhyloZooValueError."""
        # Matrix violating triangle inequality
        matrix = np.array([
            [0, 1, 5],
            [1, 0, 1],
            [5, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        co = CircularOrdering(['A', 'B', 'C'])
        with pytest.raises(PhyloZooValueError, match="must be pseudo-metric"):
            is_kalmanson(dm, co)
    
    def test_small_kalmanson_matrix(self) -> None:
        """Test Kalmanson check on small matrix."""
        # For n < 4, Kalmanson conditions are trivially satisfied
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        co = CircularOrdering(['A', 'B', 'C'])
        assert is_kalmanson(dm, co) is True
    
    def test_kalmanson_invalid_type(self) -> None:
        """Test that non-CircularOrdering circular_order raises TypeError."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        with pytest.raises(TypeError, match="must be a CircularOrdering"):
            is_kalmanson(dm, ['A', 'B'])  # type: ignore
    
    def test_kalmanson_empty_order(self) -> None:
        """Test that empty circular_order raises PhyloZooValueError."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        co = CircularOrdering([])
        with pytest.raises(PhyloZooValueError, match="cannot be empty"):
            is_kalmanson(dm, co)
    
    def test_kalmanson_very_small_matrices(self) -> None:
        """Test that matrices with < 4 elements are trivially Kalmanson."""
        # Matrix with 2 elements
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        co = CircularOrdering(['A', 'B'])
        assert is_kalmanson(dm, co) is True
        
        # Matrix with 1 element
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        co = CircularOrdering(['A'])
        assert is_kalmanson(dm, co) is True
    
    def test_kalmanson_large_matrix_performance(self) -> None:
        """Test that is_kalmanson is fast on larger matrices with numba."""
        import time
        
        # Create a large Kalmanson matrix (circular distances)
        size = 50
        matrix = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                # Circular distance
                dist = min(abs(i - j), size - abs(i - j))
                matrix[i, j] = float(dist)
        
        labels = [f'L{i}' for i in range(size)]
        dm = DistanceMatrix(matrix, labels=labels)
        circular_order = CircularOrdering(labels.copy())
        
        # Warmup
        _ = is_kalmanson(dm, circular_order)
        
        # Time the check
        start = time.perf_counter()
        result = is_kalmanson(dm, circular_order)
        elapsed = time.perf_counter() - start
        
        # Should complete in reasonable time (< 1 second for 50x50)
        assert elapsed < 1.0, f"is_kalmanson took {elapsed:.3f}s, expected < 1.0s"
        assert result is True, "Large matrix should be Kalmanson"
        
        # With numba, should be quite fast (< 0.5s)
        assert elapsed < 0.5, (
            f"is_kalmanson took {elapsed:.3f}s for {size}x{size} matrix, "
            f"expected < 0.5s with numba optimization"
        )


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

