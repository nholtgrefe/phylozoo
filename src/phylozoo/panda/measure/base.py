"""
Base functions for diversity calculations.

This module provides four common functions that work with any diversity measure:
1. diversity - compute diversity of a set of taxa
2. marginal_diversities - compute marginal diversity contributions
3. greedy_max_diversity - greedily find k taxa with highest diversity
4. solve_max_diversity - solve the maximization problem
"""

from __future__ import annotations

from typing import Any, Dict, Set

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.utils.exceptions import PhyloZooValueError, PhyloZooNotImplementedError

from .protocol import DiversityMeasure


def diversity(
    network: DirectedPhyNetwork,
    taxa: set[str],
    measure: DiversityMeasure,
    **kwargs: Any,
) -> float:
    """
    Compute the diversity of a set of taxa.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The phylogenetic network.
    taxa : set[str]
        Set of taxa to compute diversity for.
    measure : DiversityMeasure
        The diversity measure to use.
    **kwargs : Any
        Measure-specific parameters.
    
    Returns
    -------
    float
        The diversity value.
    
    Raises
    ------
    PhyloZooValueError
        If any taxa are not found in the network.
    
    Examples
    --------
    >>> from phylozoo.panda.measure import diversity, all_paths
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> net = DirectedPhyNetwork(...)
    >>> div = diversity(net, {"A", "B", "C"}, all_paths)
    """
    # Validate that all taxa are in the network
    network_taxa = set(network.taxa)
    if not taxa.issubset(network_taxa):
        missing = taxa - network_taxa
        raise PhyloZooValueError(f"Taxa not found in network: {missing}")
    
    return measure.compute_diversity(network, taxa, **kwargs)


def marginal_diversities(
    network: DirectedPhyNetwork,
    saved_taxa: set[str],
    measure: DiversityMeasure,
    **kwargs: Any,
) -> Dict[str, float]:
    """
    Compute marginal diversity contributions for all taxa.
    
    For a taxon in saved_taxa, the marginal diversity is the decrease
    when removing it from the set. For a taxon not in saved_taxa,
    the marginal diversity is the increase when adding it to the set.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The phylogenetic network.
    saved_taxa : set[str]
        Currently saved taxa.
    measure : DiversityMeasure
        The diversity measure to use.
    **kwargs : Any
        Measure-specific parameters.
    
    Returns
    -------
    Dict[str, float]
        Dictionary mapping each taxon to its marginal diversity.
    
    Examples
    --------
    >>> from phylozoo.panda.measure import marginal_diversities, all_paths
    >>> marginals = marginal_diversities(net, {"A", "B"}, all_paths)
    >>> marginals["C"]  # Marginal diversity of taxon C
    0.5
    """
    total_div = diversity(network, saved_taxa, measure, **kwargs)
    marginal: Dict[str, float] = {}
    
    all_taxa = set(network.taxa)
    
    for taxon in all_taxa:
        if taxon in saved_taxa:
            # Marginal = decrease when removing
            div_minus = diversity(
                network, saved_taxa - {taxon}, measure, **kwargs
            )
            marginal[taxon] = div_minus - total_div
        else:
            # Marginal = increase when adding
            div_plus = diversity(
                network, saved_taxa | {taxon}, measure, **kwargs
            )
            marginal[taxon] = div_plus - total_div
    
    return marginal


def greedy_max_diversity(
    network: DirectedPhyNetwork,
    k: int,
    measure: DiversityMeasure,
    **kwargs: Any,
) -> tuple[float, Set[str]]:
    """
    Greedily find k taxa with highest diversity.
    
    This function iteratively selects the taxon with the highest
    marginal diversity contribution.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The phylogenetic network.
    k : int
        Number of taxa to select.
    measure : DiversityMeasure
        The diversity measure to use.
    **kwargs : Any
        Measure-specific parameters.
    
    Returns
    -------
    tuple[float, Set[str]]
        Tuple of (diversity_value, solution_set).
    
    Raises
    ------
    PhyloZooValueError
        If k is invalid (negative or greater than number of taxa).
    
    Examples
    --------
    >>> from phylozoo.panda.measure import greedy_max_diversity, all_paths
    >>> value, solution = greedy_max_diversity(net, k=5, measure=all_paths)
    >>> len(solution)
    5
    """
    all_taxa = set(network.taxa)
    
    if k < 0:
        raise PhyloZooValueError(f"k must be non-negative, got {k}")
    if k > len(all_taxa):
        raise PhyloZooValueError(f"k must be <= {len(all_taxa)}, got {k}")
    
    saved_taxa: set[str] = set()
    
    for _ in range(k):
        marginal = marginal_diversities(
            network, saved_taxa, measure, **kwargs
        )
        
        # Filter to only unsaved taxa
        candidates = {
            t: val for t, val in marginal.items()
            if t not in saved_taxa
        }
        
        if not candidates:
            break
        
        # Select best candidate
        best_taxon = max(candidates, key=candidates.get)
        saved_taxa.add(best_taxon)
    
    # Compute final diversity value
    div_value = diversity(network, saved_taxa, measure, **kwargs)
    
    return div_value, saved_taxa


def solve_max_diversity(
    network: DirectedPhyNetwork,
    k: int,
    measure: DiversityMeasure,
    **kwargs: Any,
) -> tuple[float, Set[str]]:
    """
    Solve the maximum diversity problem.
    
    Uses measure-specific optimization if available.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The phylogenetic network.
    k : int
        Number of taxa to select.
    measure : DiversityMeasure
        The diversity measure to use.
    **kwargs : Any
        Measure-specific parameters (e.g., tree_extension for all_paths).
    
    Returns
    -------
    tuple[float, Set[str]]
        Tuple of (diversity_value, solution_set).
    
    Raises
    ------
    PhyloZooNotImplementedError
        If the measure doesn't implement solve_maximization.
    
    Examples
    --------
    >>> from phylozoo.panda.measure import solve_max_diversity, all_paths
    >>> # Optimal solution (if measure supports it)
    >>> value, solution = solve_max_diversity(
    ...     net, k=5, measure=all_paths, tree_extension="optimal_XP"
    ... )
    """
    return measure.solve_maximization(network, k=k, **kwargs)

