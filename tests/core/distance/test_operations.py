"""
Tests for distance matrix operations (TSP functions).
"""

import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.distance.operations import optimal_tsp_tour, approximate_tsp_tour
from phylozoo.core.primitives.circular_ordering import CircularOrdering


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
        
        # Check that tour is a CircularOrdering
        assert isinstance(tour, CircularOrdering)
        
        # Check that tour contains all labels
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)
        
        # Verify optimal tour: ('A', 'B', 'C', 'D') has distance 6.0
        # This is optimal because any other ordering would have at least one edge of distance 2 or 3
        expected_tour = CircularOrdering(['A', 'B', 'C', 'D'])
        assert tour == expected_tour
        
        # Verify the tour distance is optimal (6.0)
        tour_distance = sum(
            dm.get_distance(tour.order[i], tour.order[(i + 1) % len(tour)])
            for i in range(len(tour))
        )
        assert abs(tour_distance - 6.0) < 1e-10
    
    def test_single_label_tour(self) -> None:
        """Test TSP tour with single label."""
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour = optimal_tsp_tour(dm)
        assert isinstance(tour, CircularOrdering)
        assert tour.order == ('A',)
    
    def test_two_label_tour(self) -> None:
        """Test TSP tour with two labels."""
        matrix = np.array([[0, 1], [1, 0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B'])
        tour = optimal_tsp_tour(dm)
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == 2
        assert set(tour.order) == {'A', 'B'}
        
        # Verify optimal tour: ('A', 'B') has distance 2.0 (1 + 1)
        # This is the only possible tour (up to rotation/reversal)
        expected_tour = CircularOrdering(['A', 'B'])
        assert tour == expected_tour
        
        # Verify the tour distance is optimal (2.0)
        tour_distance = sum(
            dm.get_distance(tour.order[i], tour.order[(i + 1) % len(tour)])
            for i in range(len(tour))
        )
        assert abs(tour_distance - 2.0) < 1e-10
    
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
        
        # Check that tour is a CircularOrdering
        assert isinstance(tour, CircularOrdering)
        
        # All labels present exactly once
        assert len(tour) == len(set(tour.order))
        assert set(tour.order) == set(dm.labels)
        
        # Verify optimal tour: ('A', 'B', 'C', 'D') has distance 8.0
        # This is optimal: A->B (2), B->C (1), C->D (1), D->A (4) = 8.0
        expected_tour = CircularOrdering(['A', 'B', 'C', 'D'])
        assert tour == expected_tour
        
        # Verify the tour distance is optimal (8.0)
        tour_distance = sum(
            dm.get_distance(tour.order[i], tour.order[(i + 1) % len(tour)])
            for i in range(len(tour))
        )
        assert abs(tour_distance - 8.0) < 1e-10


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
        
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)
    
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
        
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)
    
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
        
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == len(dm.labels)
        assert set(tour.order) == set(dm.labels)
    
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
        assert isinstance(tour, CircularOrdering)
        assert tour.order == ('A',)
    
    def test_empty_matrix_handling(self) -> None:
        """Test handling of edge cases."""
        # Single label
        matrix = np.array([[0]])
        dm = DistanceMatrix(matrix, labels=['A'])
        tour1 = approximate_tsp_tour(dm, method='greedy')
        tour2 = optimal_tsp_tour(dm)
        assert isinstance(tour1, CircularOrdering)
        assert isinstance(tour2, CircularOrdering)
        assert tour1.order == ('A',)
        assert tour2.order == ('A',)
    
    def test_empty_matrix(self) -> None:
        """Test handling of empty distance matrix."""
        matrix = np.zeros((0, 0))
        dm = DistanceMatrix(matrix, labels=[])
        tour1 = approximate_tsp_tour(dm, method='greedy')
        tour2 = optimal_tsp_tour(dm)
        assert isinstance(tour1, CircularOrdering)
        assert isinstance(tour2, CircularOrdering)
        assert len(tour1) == 0
        assert len(tour2) == 0


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
            assert isinstance(tour, CircularOrdering)
            assert set(tour.order) == set(dm.labels)
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
        
        assert isinstance(optimal, CircularOrdering)
        assert isinstance(approximate, CircularOrdering)
        assert set(optimal.order) == set(approximate.order)
        assert set(optimal.order) == set(dm.labels)
    
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
        
        assert isinstance(tour1, CircularOrdering)
        assert isinstance(tour2, CircularOrdering)
        # CircularOrderings are equal if they represent the same circular ordering
        assert tour1 == tour2
        
        # Verify optimal tour: ('A', 'B', 'C') has distance 4.0
        # This is optimal: A->B (1), B->C (1), C->A (2) = 4.0
        # Alternative ('A', 'C', 'B') would be: A->C (2), C->B (1), B->A (1) = 4.0 (same)
        # But ('A', 'B', 'C') is the canonical form
        expected_tour = CircularOrdering(['A', 'B', 'C'])
        assert tour1 == expected_tour
        
        # Verify the tour distance is optimal (4.0)
        tour_distance = sum(
            dm.get_distance(tour1.order[i], tour1.order[(i + 1) % len(tour1)])
            for i in range(len(tour1))
        )
        assert abs(tour_distance - 4.0) < 1e-10
    
    def test_larger_matrix_approximate(self) -> None:
        """Test approximate TSP on larger matrix."""
        n = 8
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                matrix[i, j] = abs(i - j)
        dm = DistanceMatrix(matrix)
        
        tour = approximate_tsp_tour(dm, method='greedy')
        assert isinstance(tour, CircularOrdering)
        assert len(tour) == n
        assert set(tour.order) == set(range(n))

