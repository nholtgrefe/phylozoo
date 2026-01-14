"""
Unresolve tree module.

This module implements split support computation and tree unresolution by
iteratively contracting least-supported splits.
"""

import itertools
from collections.abc import Iterator
from typing import TYPE_CHECKING, Any

from ...core.network.sdnetwork import SemiDirectedPhyNetwork
from ...core.network.sdnetwork.classifications import is_tree
from ...core.network.sdnetwork.conversions import sdnetwork_from_graph
from ...core.network.sdnetwork.derivations import split_from_cutedge
from ...core.network.sdnetwork.features import cut_edges
from ...core.quartet.qprofileset import QuartetProfileSet
from ...core.split.base import Split
from ...utils.exceptions import PhyloZooValueError

if TYPE_CHECKING:
    pass


def split_support(
    profileset: QuartetProfileSet,
    split: Split,
) -> float:
    """
    Compute the split support for a given split.
    
    Given a split A|B, returns the (weighted) ratio of quartets on {a1, a2, b1, b2}
    (where a1, a2 ∈ A and b1, b2 ∈ B) that agree with the split.
    Only considers trivial, resolved profiles (single quartet with a split).
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set.
    split : Split
        The split to compute support for. Must be non-trivial and on the same
        taxa as the profile set.
    
    Returns
    -------
    float
        The split support (between 0 and 1), or 0.0 if no quartets are found.
    
    Raises
    ------
    PhyloZooValueError
        If the split is trivial or if its taxa don't match the profile set.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> split = Split({'A', 'B'}, {'C', 'D'})
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
    >>> profileset = QuartetProfileSet(profiles=[q1, q2])
    >>> support = split_support(profileset, split)
    >>> 0.0 <= support <= 1.0
    True
    """
    if split.is_trivial():
        raise PhyloZooValueError("Split must be non-trivial")
    
    if split.elements != profileset.taxa:
        raise PhyloZooValueError("Split taxa must match profile set taxa")
    
    support = 0.0
    total_weight = 0.0
    
    # Iterate over all pairs from set1 and set2
    for a1, a2 in itertools.combinations(split.set1, 2):
        for b1, b2 in itertools.combinations(split.set2, 2):
            four_taxa = frozenset({a1, a2, b1, b2})
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
            total_weight += profile_weight * quartet_weight
            
            # Check if the quartet split matches the given split
            if quartet_split == Split({a1, a2}, {b1, b2}):
                support += profile_weight * quartet_weight
    
    if total_weight == 0:
        return 0.0
    
    return support / total_weight


def unresolve_tree(
    tree: SemiDirectedPhyNetwork,
    profileset: QuartetProfileSet,
) -> Iterator[SemiDirectedPhyNetwork]:
    """
    Generate a sequence of trees by iteratively contracting least-supported splits.
    
    Starting from the given tree, computes split support for all non-trivial splits,
    sorts them by support (lowest first), and iteratively contracts them. Yields
    the tree after each contraction step.
    
    Parameters
    ----------
    tree : SemiDirectedPhyNetwork
        Starting tree. Must be a valid tree with the same taxa as the profile set.
    profileset : QuartetProfileSet
        The quartet profile set used to compute split support.
    
    Yields
    ------
    SemiDirectedPhyNetwork
        The tree after each contraction step.
    
    Raises
    ------
    PhyloZooValueError
        If the tree is not a valid tree, or if its taxa don't match the profile set.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> # Create a tree
    >>> tree = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> # Create profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> profileset = QuartetProfileSet(profiles=[q1])
    >>> # Generate contraction sequence
    >>> for contracted_tree in unresolve_tree(tree, profileset):
    ...     assert is_tree(contracted_tree)
    ...     break  # Just test first iteration
    """
    if not is_tree(tree):
        raise PhyloZooValueError("Tree must be a valid tree")
    
    if tree.taxa != profileset.taxa:
        raise PhyloZooValueError("Tree taxa must match profile set taxa")
    
    current_tree = tree
    
    # Get all cut-edges and compute their splits from the initial tree
    cut_edges_set = cut_edges(current_tree)
    split_to_edge: dict[Split, tuple[Any, Any, int]] = {}
    for u, v, key in cut_edges_set:
        edge_split = split_from_cutedge(current_tree, u, v, key=key)
        split_to_edge[edge_split] = (u, v, key)
    
    # Filter to non-trivial splits only
    non_trivial_splits = [split for split in split_to_edge.keys() if not split.is_trivial()]
    
    if not non_trivial_splits:
        # No splits to contract, yield the original tree and return
        yield current_tree
        return
    
    # Compute split support for all non-trivial splits once at the start
    split_supports: dict[Split, float] = {}
    for split in non_trivial_splits:
        split_supports[split] = split_support(profileset, split)
    
    # Sort splits by support (lowest first - least supported first)
    sorted_splits = sorted(non_trivial_splits, key=lambda x: split_supports[x])
    
    # Yield the original tree first
    yield current_tree
    
    # Iteratively contract splits in order of support
    for split in sorted_splits:      
        u, v, key = split_to_edge[split]
        
        # Contract the split by merging the endpoints (v into u)
        graph_copy = current_tree._graph.copy()
        from ...core.primitives.m_multigraph.transformations import identify_vertices
        identify_vertices(graph_copy, [u, v])
        current_tree = sdnetwork_from_graph(graph_copy, network_type='semi-directed')
        
        # Remove the contracted edge from the dictionary
        del split_to_edge[split]
        
        # Update all other edges: replace v with u (since v is merged into u)
        updated_split_to_edge: dict[Split, tuple[Any, Any, int]] = {}
        for s, (eu, ev, ek) in split_to_edge.items():
            # Replace v with u in edge endpoints
            new_u = u if eu == v else eu
            new_v = u if ev == v else ev
            # Skip self-loops (edges that would become (u, u))
            if new_u == new_v:
                continue
            updated_split_to_edge[s] = (new_u, new_v, ek)
        
        split_to_edge = updated_split_to_edge
        
        yield current_tree


