"""
Tests for distance matrix operations (TSP functions).
"""

import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.distance.operations import optimal_tsp_tour, approximate_tsp_tour


class TestOptimalTspTour:
    """Test optimal_tsp_tour function."""
    
    def test_small_tsp_tour(self) -> None:
        """Test optimal TSP tour on small matrix."""
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = optimal_tsp_tour(dm)
        
        # Check that tour contains all labels
        assert len(tour) == len(dm.labels)
        assert set(tour) == set(dm.labels)
        assert tour[0] == dm.labels[0]  # Starts at first label
    
    def test_single_label_tour(self) -> None:
        """Test TSP tour with single label."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour = optimal_tsp_tour(dm)
        assert tour == ['A']
    
    def test_two_label_tour(self) -> None:
        """Test TSP tour with two labels."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        tour = optimal_tsp_tour(dm)
        assert len(tour) == 2
        assert set(tour) == {'A', 'B'}
        assert tour[0] == 'A'  # Starts at first label
    
    def test_optimal_tour_property(self) -> None:
        """Test that optimal tour visits all labels exactly once."""
        matrix = np.array([
            [0, 2, 3, 4],
            [2, 0, 1, 2],
            [3, 1, 0, 1],
            [4, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = optimal_tsp_tour(dm)
        
        # All labels present exactly once
        assert len(tour) == len(set(tour))
        assert set(tour) == set(dm.labels)


class TestApproximateTspTour:
    """Test approximate_tsp_tour function."""
    
    def test_simulated_annealing(self) -> None:
        """Test simulated annealing method."""
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = approximate_tsp_tour(dm, method='simulated_annealing')
        
        assert len(tour) == len(dm.labels)
        assert set(tour) == set(dm.labels)
    
    def test_greedy_method(self) -> None:
        """Test greedy method."""
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = approximate_tsp_tour(dm, method='greedy')
        
        assert len(tour) == len(dm.labels)
        assert set(tour) == set(dm.labels)
    
    def test_christofides_method(self) -> None:
        """Test Christofides method."""
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        tour = approximate_tsp_tour(dm, method='christofides')
        
        assert len(tour) == len(dm.labels)
        assert set(tour) == set(dm.labels)
    
    def test_invalid_method_raises_error(self) -> None:
        """Test that invalid method raises ValueError."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        with pytest.raises(ValueError, match="Method must be one of"):
            approximate_tsp_tour(dm, method='invalid_method')
    
    def test_single_label_approximate(self) -> None:
        """Test approximate TSP with single label."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour = approximate_tsp_tour(dm, method='greedy')
        assert tour == ['A']
    
    def test_empty_matrix_handling(self) -> None:
        """Test handling of edge cases."""
        # Single label
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour1 = approximate_tsp_tour(dm, method='greedy')
        tour2 = optimal_tsp_tour(dm)
        assert tour1 == ['A']
        assert tour2 == ['A']


class TestTspTourProperties:
    """Test properties of TSP tours."""
    
    def test_tour_contains_all_labels(self) -> None:
        """Test that tour contains all labels."""
        matrix = np.array([
            [0, 1, 2, 3, 4],
            [1, 0, 1, 2, 3],
            [2, 1, 0, 1, 2],
            [3, 2, 1, 0, 1],
            [4, 3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        
        for method in ['simulated_annealing', 'greedy', 'christofides']:
            tour = approximate_tsp_tour(dm, method=method)
            assert set(tour) == set(dm.labels)
            assert len(tour) == len(dm.labels)
    
    def test_optimal_vs_approximate(self) -> None:
        """Test that optimal and approximate tours have same labels."""
        matrix = np.array([
            [0, 1, 2, 3],
            [1, 0, 1, 2],
            [2, 1, 0, 1],
            [3, 2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        
        optimal = optimal_tsp_tour(dm)
        approximate = approximate_tsp_tour(dm, method='greedy')
        
        assert set(optimal) == set(approximate)
        assert set(optimal) == set(dm.labels)
    
    def test_deterministic_optimal(self) -> None:
        """Test that optimal tour is deterministic."""
        matrix = np.array([
            [0, 1, 2],
            [1, 0, 1],
            [2, 1, 0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        
        tour1 = optimal_tsp_tour(dm)
        tour2 = optimal_tsp_tour(dm)
        
        assert tour1 == tour2
    
    def test_larger_matrix_approximate(self) -> None:
        """Test approximate TSP on larger matrix."""
        n = 8
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                matrix[i, j] = abs(i - j)
        dm = DistanceMatrix(matrix)
        
        tour = approximate_tsp_tour(dm, method='greedy')
        assert len(tour) == n
        assert set(tour) == set(range(n))

