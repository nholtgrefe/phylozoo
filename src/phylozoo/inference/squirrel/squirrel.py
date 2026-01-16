"""
Squirrel algorithm module.

This module provides the main Squirrel algorithm for reconstructing phylogenetic networks
from quartet profiles.
"""

from typing import TYPE_CHECKING, Any

from ...core.network.dnetwork import DirectedPhyNetwork
from ...core.network.sdnetwork import SemiDirectedPhyNetwork
from ...utils.parallel import ParallelConfig, ParallelBackend

from ..utils.rooting import root_at_outgroup
from .cycle_resolution import resolve_cycles
from .qjoining import adapted_quartet_joining
from .qsimilarity import sqprofileset_from_network, sqprofileset_similarity
from .sqprofile import SqQuartetProfile
from .sqprofileset import SqQuartetProfileSet
from .tstar_tree import tstar_tree
from .unresolve_tree import unresolve_tree

if TYPE_CHECKING:
    from ...core.network.dnetwork import DirectedPhyNetwork


def _process_contracted_tree(
    args: tuple[SemiDirectedPhyNetwork | str, SqQuartetProfileSet | str, str | None, dict[str, Any]],
) -> tuple[SemiDirectedPhyNetwork | str, float]:
    """
    Process a single contracted tree: resolve cycles and compute similarity.
    
    This is a helper function for parallelization. It can handle both:
    - Direct objects (for threading/sequential): (tree, profileset, outgroup, kwargs)
    - Serialized inputs (for multiprocessing): (tree_pzdot, profileset_pz, outgroup, kwargs)
    
    Parameters
    ----------
    args : tuple
        For threading/sequential:
        - contracted_tree: SemiDirectedPhyNetwork
        - profileset: SqQuartetProfileSet
        - outgroup: str | None
        - kwargs: dict[str, Any] for resolve_cycles
        
        For multiprocessing:
        - tree_pzdot: str (PhyloZoo-DOT string representation of contracted tree)
        - profileset_pz: str (PhyloZoo format string of SqQuartetProfileSet)
        - outgroup: str | None
        - kwargs: dict[str, Any] for resolve_cycles
    
    Returns
    -------
    tuple[SemiDirectedPhyNetwork | str, float]
        The network (or PhyloZoo-DOT string for multiprocessing) and its similarity score.
    """
    tree_input, profileset_input, outgroup, kwargs = args
    
    # Check if inputs are serialized (multiprocessing) or direct (threading/sequential)
    is_serialized = isinstance(tree_input, str)
    
    if is_serialized:
        # Multiprocessing: deserialize inputs
        # Import io module to ensure format is registered in worker process
        from ...core.network.sdnetwork import io as sdnetwork_io  # noqa: F401
        from ...core.network.sdnetwork import SemiDirectedPhyNetwork
        
        tree_pzdot = tree_input
        profileset_pz = profileset_input  # type: ignore
        
        # Reconstruct tree from PhyloZoo-DOT string using IOMixin
        contracted_tree = SemiDirectedPhyNetwork.from_string(tree_pzdot, format='phylozoo-dot')
        
        # Reconstruct profileset from PhyloZoo format string
        # Import here to ensure it's available in worker process
        # Also import io module to ensure format is registered
        from . import io  # noqa: F401
        from .sqprofileset import SqQuartetProfileSet
        profileset = SqQuartetProfileSet.from_string(profileset_pz, format='pz')
    else:
        # Threading/sequential: use objects directly
        contracted_tree = tree_input  # type: ignore
        profileset = profileset_input  # type: ignore
    
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
    
    if is_serialized:
        # Multiprocessing: serialize network to PhyloZoo-DOT string using IOMixin
        net_pzdot = net_new.to_string(format='phylozoo-dot')
        return (net_pzdot, score_new)
    else:
        # Threading/sequential: return network object directly
        return (net_new, score_new)


def squirrel(
    profileset: SqQuartetProfileSet,
    outgroup: str | None = None,
    parallel: ParallelConfig | None = None,
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
    parallel : ParallelConfig | None, optional
        Parallelization configuration. If None, uses sequential execution.
        By default None.
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
    >>> from phylozoo.utils.parallel import ParallelConfig, ParallelBackend
    >>> 
    >>> # Create a simple profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> profile1 = SqQuartetProfile([q1])
    >>> profileset = SqQuartetProfileSet(profiles=[profile1])
    >>> 
    >>> # Reconstruct network (sequential)
    >>> network = squirrel(profileset)
    >>> isinstance(network, SemiDirectedPhyNetwork)
    True
    >>> 
    >>> # Reconstruct network with parallelization
    >>> network = squirrel(
    ...     profileset,
    ...     parallel=ParallelConfig(
    ...         backend=ParallelBackend.MULTIPROCESSING,
    ...         n_jobs=4
    ...     )
    ... )
    >>> isinstance(network, SemiDirectedPhyNetwork)
    True
    """
    # Step 1: Get T* tree
    tstar = tstar_tree(profileset)
    
    # Step 2: Get quartet joining tree
    qj_tree = adapted_quartet_joining(profileset, starting_tree=tstar)
    
    # Step 3: Process each contracted tree (potentially in parallel)
    if parallel is None or parallel.backend == ParallelBackend.SEQUENTIAL:
        # Sequential: process trees directly as they're yielded
        networks: list[SemiDirectedPhyNetwork] = []
        scores: list[float] = []
        for contracted_tree in unresolve_tree(qj_tree, profileset):
            net_new, score_new = _process_contracted_tree(
                (contracted_tree, profileset, outgroup, kwargs)
            )
            networks.append(net_new)
            scores.append(score_new)
    elif parallel.backend == ParallelBackend.MULTIPROCESSING:
        # Multiprocessing: collect trees first, serialize, then process in parallel
        contracted_trees = list(unresolve_tree(qj_tree, profileset))
        
        # Serialize profileset to PhyloZoo format (JSON string, picklable)
        profileset_pz = profileset.to_string(format='pz')
        
        # Serialize trees to PhyloZoo-DOT strings using IOMixin
        process_args = [
            (tree.to_string(format='phylozoo-dot'), profileset_pz, outgroup, kwargs)
            for tree in contracted_trees
        ]
        
        # Process all trees using the executor
        executor = parallel.get_executor()
        try:
            results = list(executor.map(_process_contracted_tree, process_args))
        finally:
            # Clean up multiprocessing pool if needed
            if hasattr(executor, '_pool') and executor._pool is not None:
                executor._pool.close()
                executor._pool.join()
        
        # Deserialize networks from PhyloZoo-DOT strings using IOMixin
        networks = [
            SemiDirectedPhyNetwork.from_string(net_pzdot, format='phylozoo-dot')
            for net_pzdot, _ in results
        ]
        scores = [score for _, score in results]
    else:  # THREADING
        # Threading: collect trees first, then process in parallel
        contracted_trees = list(unresolve_tree(qj_tree, profileset))
        
        process_args = [
            (tree, profileset, outgroup, kwargs) for tree in contracted_trees
        ]
        
        executor = parallel.get_executor()
        try:
            results = list(executor.map(_process_contracted_tree, process_args))
        finally:
            # Clean up thread pool if needed
            if hasattr(executor, '_executor'):
                executor._executor.shutdown(wait=True)
        
        networks = [net for net, _ in results]
        scores = [score for _, score in results]
    
    # Step 6: Find best network (highest similarity score)
    best_network_index = scores.index(max(scores))
    best_network = networks[best_network_index]
    
    # Step 7: Convert to directed network if outgroup is specified
    if outgroup is not None:
        best_network = root_at_outgroup(best_network, outgroup)
    
    return best_network

