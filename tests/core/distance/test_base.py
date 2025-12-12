"""
Tests for DistanceMatrix base class.
"""

import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix


class TestDistanceMatrixCreation:
    """Test DistanceMatrix creation and validation."""
    
    def test_create_from_numpy_array(self) -> None:
        """Create DistanceMatrix from numpy array."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert len(dm) == 3
        assert dm.labels == (0, 1, 2)
    
    def test_create_with_labels(self) -> None:
        """Create DistanceMatrix with custom labels."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        assert len(dm) == 2
        assert dm.labels == ('A', 'B')
    
    def test_non_square_matrix_raises_error(self) -> None:
        """Non-square matrix should raise ValueError."""
        matrix = np.array([[0, 1, 2], [1, 0, 1]])
        with pytest.raises(ValueError, match="must be square"):
            DistanceMatrix(matrix)
    
    def test_non_symmetric_matrix_raises_error(self) -> None:
        """Non-symmetric matrix should raise ValueError."""
        matrix = np.array([[0, 1, 2], [2, 0, 1], [1, 2, 0]])
        with pytest.raises(ValueError, match="must be symmetric"):
            DistanceMatrix(matrix)
    
    def test_wrong_type_raises_error(self) -> None:
        """Non-numpy array should raise TypeError."""
        with pytest.raises(TypeError, match="must be a numpy ndarray"):
            DistanceMatrix([[0, 1], [1, 0]])  # type: ignore
    
    def test_labels_size_mismatch_raises_error(self) -> None:
        """Labels list size mismatch should raise ValueError."""
        matrix = np.array([[0, 1], [1, 0]])
        with pytest.raises(ValueError, match="must match"):
            DistanceMatrix(matrix, labels=['A', 'B', 'C'])
    
    def test_2d_array_required(self) -> None:
        """1D array should raise ValueError."""
        matrix = np.array([0, 1, 2])
        with pytest.raises(ValueError, match="must be a 2D array"):
            DistanceMatrix(matrix)
    
    def test_symmetric_with_tolerance(self) -> None:
        """Matrix symmetric within floating point tolerance should work."""
        matrix = np.array([[0, 1.0 + 1e-12], [1.0, 0]])
        dm = DistanceMatrix(matrix)
        assert len(dm) == 2


class TestDistanceMatrixProperties:
    """Test DistanceMatrix properties and accessors."""
    
    def test_len_property(self) -> None:
        """Test __len__ method."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix)
        assert len(dm) == 2
    
    def test_labels_property(self) -> None:
        """Test labels property."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        assert dm.labels == ('A', 'B')
        assert isinstance(dm.labels, tuple)  # Immutable
    
    def test_indices_property(self) -> None:
        """Test indices property."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert dm.indices == (0, 1, 2)
        assert isinstance(dm.indices, tuple)  # Immutable
    
    def test_np_array_property(self) -> None:
        """Test np_array property."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        arr = dm.np_array
        assert isinstance(arr, np.ndarray)
        assert arr.shape == (2, 2)
        assert arr[0, 1] == 1.0
        # Test immutability
        try:
            arr[0, 1] = 999
            assert False, "Array should be read-only"
        except ValueError:
            pass  # Expected
    
    def test_get_distance(self) -> None:
        """Test get_distance method."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        assert dm.get_distance('A', 'B') == 1.0
        assert dm.get_distance('B', 'A') == 1.0  # Symmetric
        assert dm.get_distance('A', 'C') == 2.0
    
    def test_get_distance_invalid_label(self) -> None:
        """get_distance with invalid label should raise ValueError."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        with pytest.raises(ValueError, match="not found"):
            dm.get_distance('C', 'A')
    
    def test_get_index(self) -> None:
        """Test get_index method."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        assert dm.get_index('A') == 0
        assert dm.get_index('B') == 1
    
    def test_get_index_invalid_label(self) -> None:
        """get_index with invalid label should raise ValueError."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        with pytest.raises(ValueError, match="not found"):
            dm.get_index('C')


class TestDistanceMatrixImmutability:
    """Test that DistanceMatrix is immutable."""
    
    def test_labels_are_immutable(self) -> None:
        """Elements tuple should be immutable."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        assert isinstance(dm.labels, tuple)
        # Cannot modify tuple (would raise TypeError)
        with pytest.raises((AttributeError, TypeError)):
            dm.labels.append('C')  # type: ignore
    
    def test_copy_creates_new_instance(self) -> None:
        """copy() should create a new instance."""
        matrix = np.array([[0, 1], [1, 0]])
        dm1 = DistanceMatrix(matrix, labels=['A', 'B'])
        dm2 = dm1.copy()
        assert dm1 is not dm2
        assert dm1.labels == dm2.labels
        assert dm1.get_distance('A', 'B') == dm2.get_distance('A', 'B')


class TestDistanceMatrixMagicMethods:
    """Test magic methods."""
    
    def test_len(self) -> None:
        """Test __len__ method."""
        matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        dm = DistanceMatrix(matrix)
        assert len(dm) == 3
    
    def test_contains(self) -> None:
        """Test __contains__ method."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        assert 'A' in dm
        assert 'B' in dm
        assert 'C' not in dm
    
    def test_repr(self) -> None:
        """Test __repr__ method."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        repr_str = repr(dm)
        assert 'DistanceMatrix' in repr_str
        assert 'size=2' in repr_str or 'labels=' in repr_str
    
    def test_str(self) -> None:
        """Test __str__ method."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        str_repr = str(dm)
        assert 'DistanceMatrix' in str_repr
        assert 'A' in str_repr or 'B' in str_repr


class TestDistanceMatrixEdgeCases:
    """Test edge cases."""
    
    def test_single_label(self) -> None:
        """Test DistanceMatrix with single label."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        assert len(dm) == 1
        assert dm.get_distance('A', 'A') == 0.0
    
    def test_empty_matrix_not_allowed(self) -> None:
        """Empty matrix should raise ValueError."""
        matrix = np.array([[]])
        with pytest.raises(ValueError):
            DistanceMatrix(matrix)
    
    def test_large_matrix(self) -> None:
        """Test DistanceMatrix with larger size."""
        n = 10
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                matrix[i, j] = abs(i - j)
        dm = DistanceMatrix(matrix)
        assert len(dm) == n
        assert dm.get_distance(0, n-1) == float(n-1)

