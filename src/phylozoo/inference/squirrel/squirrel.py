"""
Squirrel algorithm module.

This module provides the main Squirrel algorithm for reconstructing phylogenetic networks
from quartet profiles.
"""

from typing import TYPE_CHECKING, Any

from ...core.network.dnetwork import DirectedPhyNetwork
from ...core.network.sdnetwork import SemiDirectedPhyNetwork

from ..utils.rooting import root_at_outgroup
from .cycle_resolution import resolve_cycles
from .qjoining import adapted_quartet_joining
from .qsimilarity import sqprofileset_from_network, sqprofileset_similarity
from .sqprofileset import SqQuartetProfileSet
from .tstar_tree import tstar_tree
from .unresolve_tree import unresolve_tree

if TYPE_CHECKING:
    from ...core.network.dnetwork import DirectedPhyNetwork


def squirrel(
    profileset: SqQuartetProfileSet,
    outgroup: str | None = None,
    **kwargs: Any,
) -> SemiDirectedPhyNetwork | DirectedPhyNetwork:
    """
    Reconstruct a phylogenetic network from a squirrel quartet profile set.
    
    This is the main Squirrel algorithm. It:
    1. Computes the T* tree from the profile set
    2. Applies quartet joining to get the quartet-joining tree
    3. Uses unresolve_tree to iteratively contract least-supported splits
    4. For each contracted tree, reconstructs a network and computes similarity
    5. Returns the network with the highest similarity score
    
    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The squirrel quartet profile set. Must be dense.
    outgroup : str | None, optional
        Optional outgroup taxon. If specified, the function returns a DirectedPhyNetwork
        rooted at the edge leading to the outgroup. By default None.
    **kwargs : Any
        Additional keyword arguments passed to `resolve_cycles`:
        - rho : tuple[float, float, float, float], optional
            Rho vector (rho_c, rho_s, rho_a, rho_o) for distance computation in cycle resolution.
            By default (0.5, 1.0, 0.5, 1.0) (Squirrel/MONAD).
        - tsp_threshold : int | None, optional
            Maximum partition size for which to solve TSP optimally in cycle resolution.
            If partition size is larger, use the approximate method. If None, always use optimal.
            By default 13.
    
    Returns
    -------
    SemiDirectedPhyNetwork | DirectedPhyNetwork
        The network with the highest similarity score.
        If outgroup is None: Returns SemiDirectedPhyNetwork
        If outgroup is specified: Returns DirectedPhyNetwork
    
    Raises
    ------
    PhyloZooValueError
        If profileset is not dense.
        If outgroup is specified but not found in the network.
    
    Examples
    --------
    >>> from phylozoo.inference.squirrel.sqprofile import SqQuartetProfile
    >>> from phylozoo.inference.squirrel.sqprofileset import SqQuartetProfileSet
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.split.base import Split
    >>> 
    >>> # Create a simple profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> profile1 = SqQuartetProfile([q1])
    >>> profileset = SqQuartetProfileSet(profiles=[profile1])
    >>> 
    >>> # Reconstruct network
    >>> network = squirrel(profileset)
    >>> isinstance(network, SemiDirectedPhyNetwork)
    True
    """
    # Step 1: Get T* tree
    tstar = tstar_tree(profileset)
    
    # Step 2: Get quartet joining tree
    qj_tree = adapted_quartet_joining(profileset, starting_tree=tstar)
    
    # Step 3: Use unresolve_tree to get contraction sequence
    # This already handles split support computation and sorting
    networks: list[SemiDirectedPhyNetwork] = []
    scores: list[float] = []
    
    # Process each tree in the contraction sequence
    for contracted_tree in unresolve_tree(qj_tree, profileset):
        # Reconstruct network from the contracted tree
        net_new = resolve_cycles(
            profileset,
            contracted_tree,
            outgroup=outgroup,
            **kwargs,
        )
        
        # Compute similarity score
        quarnets_new = sqprofileset_from_network(net_new)
        score_new = sqprofileset_similarity(profileset, quarnets_new)
        
        networks.append(net_new)
        scores.append(score_new)
    
    # Step 4: Find best network (highest similarity score)
    best_network_index = scores.index(max(scores))
    best_network = networks[best_network_index]
    
    # Step 5: Convert to directed network if outgroup is specified
    if outgroup is not None:
        best_network = root_at_outgroup(best_network, outgroup)
    
    return best_network

