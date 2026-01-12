"""
Protocol defining the interface for diversity measures.

This module provides the DiversityMeasure protocol that all diversity
measure implementations must follow.
"""

from __future__ import annotations

from typing import Any, Protocol, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork


class DiversityMeasure(Protocol):
    """
    Protocol defining the interface for diversity measures.
    
    Measures can implement only the methods they need. The base functions
    will handle common algorithms using the measure's compute_diversity method.
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute the diversity of a set of taxa.
        
        Parameters
        ----------
        network : DirectedPhyNetwork
            The phylogenetic network.
        taxa : Set[str]
            Set of taxa to compute diversity for.
        **kwargs : Any
            Measure-specific parameters.
        
        Returns
        -------
        float
            The diversity value.
        
        Raises
        ------
        ValueError
            If any taxa are not found in the network.
        """
        ...

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        k: int,
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve the maximum diversity problem.
        
        Measures can implement this for measure-specific optimization
        (e.g., scanwidth-based DP for all_paths).
        
        Parameters
        ----------
        network : DirectedPhyNetwork
            The phylogenetic network.
        k : int
            Number of taxa to select.
        **kwargs : Any
            Measure-specific parameters (e.g., tree_extension for all_paths).
        
        Returns
        -------
        tuple[float, Set[str]]
            Tuple of (diversity_value, solution_set).
        
        Raises
        ------
        NotImplementedError
            If the measure doesn't support optimization.
        """
        raise NotImplementedError(
            "This measure does not implement custom optimization."
        )

