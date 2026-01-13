"""
Quartet joining algorithm module.

This module implements the adapted quartet joining algorithm for constructing
phylogenetic trees from quartet profiles.

TODO: omegbar also for non-trivial profiles.
"""

import itertools
from typing import TYPE_CHECKING, Any

from ...core.network.sdnetwork import SemiDirectedPhyNetwork
from ...core.network.sdnetwork.classifications import is_tree
from ...core.network.sdnetwork.conversions import sdnetwork_from_graph
from ...core.network.sdnetwork.derivations import split_from_cutedge
from ...core.quartet.base import Quartet
from ...core.quartet.qprofileset import QuartetProfileSet
from ...core.split.base import Split

if TYPE_CHECKING:
    pass


def _omega_bar(
    profileset: QuartetProfileSet,
    X1: frozenset[str],
    X2: frozenset[str],
    X3: frozenset[str],
    X4: frozenset[str],
) -> float:
    """
    Compute omega_bar score for quartet joining.
    
    This function computes the proportion of quartets that support the split
    X1 ∪ X2 | X3 ∪ X4, where we consider all combinations of one taxon from each set.
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set.
    X1 : frozenset[str]
        First set of taxa.
    X2 : frozenset[str]
        Second set of taxa.
    X3 : frozenset[str]
        Third set of taxa.
    X4 : frozenset[str]
        Fourth set of taxa.
    
    Returns
    -------
    float
        The omega_bar score (between 0 and 1), or 1.0 if no quartets are found.
    
    Notes
    -----
    Only considers trivial, resolved profiles (single quartet with a split).
    """
    res = 0.0
    total_sum = 0.0
    
    for (x1, x2, x3, x4) in itertools.product(X1, X2, X3, X4):
        four_taxa = frozenset({x1, x2, x3, x4})
        profile = profileset.get_profile(four_taxa)
        
        if profile is None:
            continue
        
        # Only consider trivial, resolved profiles
        if not profile.is_trivial():
            continue
        
        # Get the single quartet
        quartet = next(iter(profile.quartets))
        
        # Check if quartet is resolved (has a split)
        quartet_split = quartet.split
        if quartet_split is None:
            continue
        
        # Get profile weight
        profile_weight = profileset.get_profile_weight(four_taxa) or 1.0
        quartet_weight = profile.get_weight(quartet) or 1.0
        total_weight = profile_weight * quartet_weight
        
        total_sum += total_weight
        
        # Check if the quartet split matches X1 ∪ X2 | X3 ∪ X4
        if (quartet_split.set1 == {x1, x2} and quartet_split.set2 == {x3, x4}) or \
           (quartet_split.set2 == {x1, x2} and quartet_split.set1 == {x3, x4}):
            res += total_weight
    
    if total_sum == 0:
        return 1.0
    
    return res / total_sum


def adapted_quartet_joining(
    profileset: QuartetProfileSet,
    starting_tree: SemiDirectedPhyNetwork,
) -> SemiDirectedPhyNetwork:
    """
    Apply the adapted quartet joining algorithm to construct a tree.
    
    This is an adapted version of the quartet joining algorithm that takes an
    arbitrary starting tree (not necessarily a star tree). The algorithm iteratively
    refines the tree by selecting the best pair of neighbors to join at each center
    node, based on quartet support scores. The algorithm continues until the tree
    is binary (all internal nodes have degree 3).
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set. Only trivial, resolved profiles are considered.
    starting_tree : SemiDirectedPhyNetwork
        Starting tree. The tree must be valid and have the same taxa as the profile set.
        Unlike the original quartet joining algorithm, this can be any tree, not
        necessarily a star tree.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        The constructed tree.
    
    Raises
    ------
    ValueError
        If the starting tree is not a tree, or if its taxa don't match the profile set.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> # Create trivial quartet profiles (each profile contains a single quartet)
    >>> # For 5 taxa: A, B, C, D, E
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
    >>> q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
    >>> profileset = QuartetProfileSet(profiles=[q1, q2, q3])
    >>> # Create a starting star tree with all 5 taxa connected to a center node
    >>> center = '_center'
    >>> starting = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(center, 'A'), (center, 'B'), (center, 'C'), (center, 'D'), (center, 'E')],
    ...     nodes=[('A', {'label': 'A'}), ('B', {'label': 'B'}), ('C', {'label': 'C'}), 
    ...            ('D', {'label': 'D'}), ('E', {'label': 'E'})]
    ... )
    >>> tree = adapted_quartet_joining(profileset, starting)
    >>> is_tree(tree)
    True
    """
    # Validate starting tree
    if not is_tree(starting_tree):
        raise ValueError("Starting tree must be a valid tree")
    if starting_tree.taxa != profileset.taxa:
        raise ValueError("Starting tree taxa must match profile set taxa")
    
    current_tree = starting_tree
    
    # Calculate number of iterations needed
    # We need n - 3 - (number of non-trivial splits) iterations
    # For a star tree, number of non-trivial splits is 0
    n = len(profileset.taxa)
    if n < 4:
        # Too few taxa for quartet joining
        return current_tree
    
    # Count existing non-trivial splits
    # In a tree: non-trivial splits = total_edges - n
    total_edges = current_tree.number_of_edges()
    nr_iterations = n - 3 - (total_edges - n)
    
    # Main iteration loop
    for _ in range(nr_iterations):
        # Find center nodes (internal nodes with degree > 3)
        center_nodes = [
            node for node in current_tree.internal_nodes
            if current_tree.degree(node) > 3
        ]
        
        if not center_nodes:
            # Tree is already binary
            break
        
        # Compute scores for all pairs of neighbors at each center node
        score: dict[tuple[Any, Any, Any], float] = {}
        
        for center_node in center_nodes:
            if current_tree.degree(center_node) == 3:
                continue
            
            # Get neighbors and their induced splits
            neighbors_list = list(current_tree.neighbors(center_node))
            C: dict[Any, frozenset[str]] = {}
            
            for u in neighbors_list:
                # Get split from edge (u, center_node) and taxa on the side of u
                _, (node_u, taxa_u), _ = split_from_cutedge(
                    current_tree, u, center_node, return_node_taxa=True
                )
                # Store the taxa on the side of u
                C[u] = taxa_u
            
            # Score all pairs of neighbors
            for (u1, u2) in itertools.combinations(neighbors_list, 2):
                score[(u1, u2, center_node)] = 0.0
                
                # Consider all pairs of components (excluding u1 and u2)
                for (A1, A2) in itertools.combinations(C.values(), 2):
                    if C[u1] not in [A1, A2] and C[u2] not in [A1, A2]:
                        score[(u1, u2, center_node)] += _omega_bar(
                            profileset,
                            C[u1],
                            C[u2],
                            A1,
                            A2,
                        )
        
        if not score:
            # No valid scores found
            break
        
        # Find the best pair
        u1star, u2star, cstar = max(score, key=score.get)
        
        # Create the split by modifying the tree
        # Remove edges (cstar, u1star) and (cstar, u2star)
        tree_graph = current_tree._graph.copy()
        
        # Remove edges (trees only have undirected edges, no parallel edges)
        if tree_graph.has_edge(cstar, u1star):
            tree_graph.remove_edge(cstar, u1star)
        
        if tree_graph.has_edge(cstar, u2star):
            tree_graph.remove_edge(cstar, u2star)
        
        # Create new internal node using generate_node_ids
        w = next(tree_graph.generate_node_ids(1))
        tree_graph.add_node(w)
        
        # Add edges: (w, u1star), (w, u2star), (w, cstar)
        tree_graph.add_undirected_edge(w, u1star)
        tree_graph.add_undirected_edge(w, u2star)
        tree_graph.add_undirected_edge(w, cstar)
        
        # Convert back to SemiDirectedPhyNetwork
        current_tree = sdnetwork_from_graph(tree_graph, network_type='semi-directed')
    
    return current_tree


def quartet_joining(
    profileset: QuartetProfileSet,
) -> SemiDirectedPhyNetwork:
    """
    Apply the quartet joining algorithm to construct a tree.
    
    This is the original quartet joining algorithm that starts with a star tree
    (all taxa connected to a single center node). The algorithm iteratively refines
    the star tree by selecting the best pair of neighbors to join at each center
    node, based on quartet support scores. The algorithm continues until the tree
    is binary (all internal nodes have degree 3).
    
    This function creates a star tree from the profile set's taxa and calls
    `adapted_quartet_joining` with it as the starting tree.
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set. Only trivial, resolved profiles are considered.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        The constructed tree.
    
    Raises
    ------
    ValueError
        If there are fewer than 4 taxa in the profile set.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> # Create trivial quartet profiles for 5 taxa
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
    >>> q3 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
    >>> profileset = QuartetProfileSet(profiles=[q1, q2, q3])
    >>> tree = quartet_joining(profileset)
    >>> is_tree(tree)
    True
    >>> sorted(tree.taxa)
    ['A', 'B', 'C', 'D', 'E']
    """
    if len(profileset.taxa) < 4:
        raise ValueError("Quartet joining requires at least 4 taxa")
    
    # Create a star tree with all taxa connected to a center node
    center = "_center"
    taxa_list = sorted(profileset.taxa)
    undirected_edges = [(center, taxon) for taxon in taxa_list]
    nodes = [(taxon, {'label': taxon}) for taxon in taxa_list]
    
    starting_tree = SemiDirectedPhyNetwork(
        undirected_edges=undirected_edges,
        nodes=nodes
    )
    
    return adapted_quartet_joining(profileset, starting_tree)

