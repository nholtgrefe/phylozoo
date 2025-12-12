"""
Tests for distance matrix classification API methods.
"""

import time

import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.distance.classifications import (
    has_zero_diagonal,
    is_nonnegative,
    satisfies_triangle_inequality,
)


class TestHasZeroDiagonal:
    """Test has_zero_diagonal function."""
    
    def test_zero_diagonal(self) -> None:
        """Test matrix with zero diagonal."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert has_zero_diagonal(dm) is True
    
    def test_non_zero_diagonal(self) -> None:
        """Test matrix with non-zero diagonal."""
        matrix = np.array([[1, 1], [1, 1]])
        dm = DistanceMatrix(matrix)
        assert has_zero_diagonal(dm) is False
    
    def test_single_element_zero(self) -> None:
        """Test single element matrix with zero."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix)
        assert has_zero_diagonal(dm) is True


class TestIsNonnegative:
    """Test is_nonnegative function."""
    
    def test_nonnegative_matrix(self) -> None:
        """Test matrix with all non-negative values."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert is_nonnegative(dm) is True
    
    def test_negative_values(self) -> None:
        """Test matrix with negative values."""
        # Note: DistanceMatrix does not validate non-negativity, only symmetry
        # So we can create a matrix with negative values
        matrix = np.array([[0, -1], [-1, 0]])
        dm = DistanceMatrix(matrix)
        # is_nonnegative should return False for this matrix
        assert is_nonnegative(dm) is False
    
    def test_zero_matrix(self) -> None:
        """Test zero matrix."""
        matrix = np.zeros((3, 3))
        dm = DistanceMatrix(matrix)
        assert is_nonnegative(dm) is True


class TestSatisfiesTriangleInequality:
    """Test satisfies_triangle_inequality function."""
    
    def test_satisfies_triangle_inequality(self) -> None:
        """Test matrix that satisfies triangle inequality."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert satisfies_triangle_inequality(dm) is True
    
    def test_violates_triangle_inequality(self) -> None:
        """Test matrix that violates triangle inequality."""
        # d(A,C) = 5 > d(A,B) + d(B,C) = 1 + 1 = 2
        matrix = np.array([[0, 1, 5], [1, 0, 1], [5, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert satisfies_triangle_inequality(dm) is False
    
    def test_single_element(self) -> None:
        """Test single element matrix (trivially satisfies)."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix)
        assert satisfies_triangle_inequality(dm) is True
    
    def test_zero_matrix(self) -> None:
        """Test zero matrix (satisfies triangle inequality)."""
        matrix = np.zeros((3, 3))
        dm = DistanceMatrix(matrix)
        assert satisfies_triangle_inequality(dm) is True
    
    def test_large_matrix_performance(self) -> None:
        """Test that triangle inequality check is fast on large matrices."""
        # Create a large valid distance matrix (Euclidean distances from random points)
        # This ensures it satisfies triangle inequality
        size = 200
        np.random.seed(42)
        points = np.random.rand(size, 2) * 10
        
        # Use vectorized computation for faster matrix creation
        matrix = np.zeros((size, size))
        for i in range(size):
            for j in range(i, size):
                dist = np.linalg.norm(points[i] - points[j])
                matrix[i, j] = dist
                matrix[j, i] = dist  # Symmetric
        
        dm = DistanceMatrix(matrix)
        
        # Warmup: first call compiles numba function
        _ = satisfies_triangle_inequality(dm)
        
        # Time the check (should be fast due to numba optimization)
        start = time.perf_counter()
        result = satisfies_triangle_inequality(dm)
        elapsed = time.perf_counter() - start
        
        # Should complete in reasonable time (e.g., < 1 second for 200x200)
        # With numba, this should be much faster (< 0.1s typically)
        assert elapsed < 1.0, f"Triangle inequality check took {elapsed:.3f}s, expected < 1.0s"
        assert result is True, "Large matrix should satisfy triangle inequality"
        
        # For a 200x200 matrix with numba, we expect it to be quite fast
        # Set a reasonable limit (0.5s) to ensure numba is working
        assert elapsed < 0.5, (
            f"Triangle inequality check took {elapsed:.3f}s for {size}x{size} matrix, "
            f"expected < 0.5s with numba optimization"
        )

