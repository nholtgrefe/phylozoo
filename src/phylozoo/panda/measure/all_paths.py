"""
All-paths diversity measure implementation.

This module implements the all-paths diversity measure, which computes
the sum of branch lengths (edge weights) in the subnetwork induced by
a set of taxa.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from typing import Any, Set

import math

import networkx as nx

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import has_parallel_edges, is_binary
from phylozoo.core.network.dnetwork.features import k_blobs
from phylozoo.core.network.dnetwork.transformations import binary_resolution
from phylozoo.panda.utils.scanwidth import ScanwidthDAG, Extension, TreeExtension
from phylozoo.utils.exceptions import PhyloZooValueError, PhyloZooNotImplementedError, PhyloZooRuntimeError

from .protocol import DiversityMeasure


def _powerset(s: set[Any]) -> Any:
    """
    Yield subsets of a set s as sets, one at a time.
    
    Parameters
    ----------
    s : set[Any]
        Set to generate powerset for.
    
    Yields
    ------
    set[Any]
        Each subset of s.
    """
    s_list = tuple(s)
    for r in range(len(s_list) + 1):
        for subset in combinations(s_list, r):
            yield set(subset)


class _DPInstance:
    """
    Internal class for running the scanwidth FPT-algorithm for MAPPD.
    
    Takes as input a binary network, a parameter k, and a tree extension.
    Implements dynamic programming over the tree extension to solve MAPPD.
    
    Attributes
    ----------
    network : DirectedPhyNetwork
        The binary phylogenetic network.
    tree_extension : TreeExtension
        The tree extension for dynamic programming.
    minus_infinity : float
        A value smaller than any possible diversity (used as sentinel).
    GW : dict[Any, set[tuple[Any, Any]]]
        Table storing sets of edges in each scanwidth bag GW_v.
    table : dict[Any, dict[int, dict[frozenset[tuple[Any, Any]], tuple[float, Any]]]]
        Dynamic programming table: [vertex][l][phi] -> (diversity_score, pointer).
    m : int
        Maximum edge offspring count (for heuristic pruning).
    """
    
    def __init__(self, network: DirectedPhyNetwork, tree_extension: TreeExtension) -> None:
        """
        Initialize DP instance.
        
        Parameters
        ----------
        network : DirectedPhyNetwork
            The binary phylogenetic network.
        tree_extension : TreeExtension
            The tree extension for dynamic programming.
        """
        self.network = network
        self.tree_extension = tree_extension
        
        # Compute minus_infinity: smaller than any possible diversity
        total_weight = 0.0
        for u, v, key, data in network._graph.edges(keys=True, data=True):
            branch_length = data.get('branch_length')
            if branch_length is not None:
                total_weight += branch_length
            else:
                total_weight += 1.0
        self.minus_infinity = -2 * total_weight - 1
        
        self.GW = self._initialize_GW()
        self.table: dict[Any, dict[int, dict[frozenset[tuple[Any, Any]], tuple[float, Any]]]] = (
            defaultdict(lambda: defaultdict(dict))
        )
        self.m = self._max_edge_offspring_count()
    
    def _get_branch_length(self, u: Any, v: Any) -> float:
        """
        Get branch length of edge (u, v), defaulting to 1.0 if not specified.
        
        Parameters
        ----------
        u : Any
            Source node.
        v : Any
            Target node.
        
        Returns
        -------
        float
            Branch length (or 1.0 if not specified).
        """
        branch_length = self.network.get_branch_length(u, v)
        if branch_length is not None:
            return branch_length
        return 1.0
    
    def _initialize_GW(self) -> dict[Any, set[tuple[Any, Any]]]:
        """
        Compute all GW_v scanwidth bags.
        
        Returns
        -------
        dict[Any, set[tuple[Any, Any]]]
            Dictionary mapping each node to its scanwidth bag (set of edges).
        """
        res: dict[Any, set[tuple[Any, Any]]] = {v: set() for v in self.network._graph.nodes()}
        
        for v in self.network._graph.nodes():
            sink = set(nx.descendants(self.tree_extension.tree, v))
            sink.add(v)
            
            for u, w, key in self.network._graph.edges(keys=True):
                if u not in sink and w in sink:
                    res[v].add((u, w))
        
        return res
    
    def _max_edge_offspring_count(self) -> int:
        """
        Return the maximum value m_x over all leaves x of the network,
        where m_x = |{e: x \in off(e)}|.
        
        Returns
        -------
        int
            Maximum edge offspring count.
        """
        # Indegree map
        indeg = {v: self.network.indegree(v) for v in self.network._graph.nodes()}
        
        # Memoization
        edge_count: dict[Any, int] = {}
        
        def count_edges(u: Any) -> int:
            """Recursively count edges in subnetwork below u."""
            if u in edge_count:
                return edge_count[u]
            total = indeg[u]
            for v in self.network.parents(u):
                total += count_edges(v)
            edge_count[u] = total
            return total
        
        best_count = 0
        for leaf in self.network.leaves:
            c = count_edges(leaf)
            if c > best_count:
                best_count = c
        
        return best_count
    
    def _fill_dp_table(self, k: int) -> None:
        """
        Bottom-up dynamic programming that fills self.table.
        
        Parameters
        ----------
        k : int
            Number of taxa to select.
        """
        self.table = defaultdict(lambda: defaultdict(dict))
        
        for v in nx.dfs_postorder_nodes(self.tree_extension.tree):
            if v in self.network.leaves:
                self._process_leaf_node(v, k)
            else:
                self._process_non_leaf_node(v, k)
    
    def _process_leaf_node(self, v: Any, k: int) -> None:
        """
        Process a leaf node v and write values to the table.
        
        Parameters
        ----------
        v : Any
            Leaf node to process.
        k : int
            Maximum number of taxa to select.
        """
        GW_v = self.GW.get(v, set())
        
        for phi in _powerset(GW_v):
            phi_frozen = frozenset(phi)
            for l in range(k + 1):
                if l == 0:
                    dp = 0.0 if len(phi) == 0 else self.minus_infinity
                    pointer = None
                else:
                    dp = sum(self._get_branch_length(u, v) for (u, v) in phi)
                    pointer = (0,)
                
                self.table[v][l][phi_frozen] = (dp, pointer)
    
    def _process_non_leaf_node(self, v: Any, k: int) -> None:
        """
        Process a non-leaf node v and write values to the table.
        
        Parameters
        ----------
        v : Any
            Non-leaf node to process.
        k : int
            Maximum number of taxa to select.
        """
        GW_v = self.GW.get(v, set())
        
        # delta_in and delta_out
        delta_in_v = set()
        for u, _, _ in self.network.incident_parent_edges(v, keys=True):
            delta_in_v.add((u, v))
        
        delta_out_v = set()
        for _, w, _ in self.network.incident_child_edges(v, keys=True):
            delta_out_v.add((v, w))
        
        # Children of v in the tree-extension
        children = list(self.tree_extension.tree.successors(v))
        GW_children = [self.GW.get(child, set()) for child in children]
        
        for phi in _powerset(GW_v):
            phi_frozen = frozenset(phi)
            phi_len = len(phi)
            
            # Used for heuristic pruning
            bound = math.ceil(phi_len / self.m) if self.m > 0 else 0
            
            # val3 equals Omega(phi \cap \delta_in (v))
            val3 = sum(
                self._get_branch_length(u, v) for (u, v) in phi & delta_in_v
            )
            
            # phi_psi_subsets contains all combinations of phi \cup psi with psi \subseteq \delta_out (v)
            phi_psi_subsets = [phi | psi for psi in _powerset(delta_out_v)]
            
            # Create list containing one list for each child x of v, containing
            # all combinations of (phi \cup psi) \cap GW_x
            phi_psi_GW_subsets: list[list[frozenset[tuple[Any, Any]]]] = []
            for i in range(len(children)):
                GW = GW_children[i]
                phi_psi_GW_subsets.append([])
                for phi_psi in phi_psi_subsets:
                    phi_psi_GW_subsets[i].append(frozenset(phi_psi & GW))
            
            for l in range(k + 1):
                # Base case: l = 0
                if l == 0:
                    dp = 0.0 if phi_len == 0 else self.minus_infinity
                    pointer = None
                
                # Heuristic pruning: if l < bound, we can never have |off(e) \cap A| > 0 for all e \in phi
                elif l < bound:
                    dp = self.minus_infinity
                    pointer = None
                
                else:
                    dp = self.minus_infinity - 1
                    pointer = None
                    
                    # Case: v has one child
                    if len(children) == 1:
                        u = children[0]
                        
                        for S in phi_psi_GW_subsets[0]:
                            val1 = self.table.get(u, {}).get(l, {}).get(S, (self.minus_infinity, None))[0]
                            val = val1 + val3
                            if val >= dp:
                                dp = val
                                pointer = (1, (u, l, S))
                    
                    # Case: v has two children
                    elif len(children) == 2:
                        u, w = children
                        
                        for l_prime in range(l + 1):
                            for i in range(len(phi_psi_GW_subsets[0])):
                                S1 = phi_psi_GW_subsets[0][i]
                                val1 = self.table.get(u, {}).get(l_prime, {}).get(S1, (self.minus_infinity, None))[0]
                                
                                S2 = phi_psi_GW_subsets[1][i]
                                val2 = self.table.get(w, {}).get(l - l_prime, {}).get(S2, (self.minus_infinity, None))[0]
                                
                                val = val1 + val2 + val3
                                if val >= dp:
                                    dp = val
                                    pointer = (2, (u, l_prime, S1), (w, l - l_prime, S2))
                
                self.table[v][l][phi_frozen] = (dp, pointer)
    
    def _backtrack_solution(self, v: Any, l: int, phi_frozen: frozenset[tuple[Any, Any]]) -> set[Any]:
        """
        Use backtracking to get the solution set corresponding to DP[v][l][phi_frozen].
        
        Parameters
        ----------
        v : Any
            Node in tree extension.
        l : int
            Number of taxa to select.
        phi_frozen : frozenset[tuple[Any, Any]]
            Frozen set of edges.
        
        Returns
        -------
        set[Any]
            Set of leaf node IDs in the solution.
        """
        entry = self.table.get(v, {}).get(l, {}).get(phi_frozen, (None, None))
        pointer = entry[1]
        
        # No saving / invalid options -> return empty set
        if pointer is None:
            return set()
        
        # v is a leaf -> return {v}
        if pointer[0] == 0:
            return {v}
        
        # v has one child -> return solutions of that child
        elif pointer[0] == 1:
            child, l_child, phi_child = pointer[1]
            return self._backtrack_solution(child, l_child, phi_child)
        
        # v has two children -> return union of solutions
        elif pointer[0] == 2:
            (u, l1, phi_u) = pointer[1]
            (w, l2, phi_w) = pointer[2]
            sol1 = self._backtrack_solution(u, l1, phi_u)
            sol2 = self._backtrack_solution(w, l2, phi_w)
            return sol1 | sol2
    
    def solve(self, k: int) -> tuple[float, set[Any]]:
        """
        Solve MAPPD with parameter k.
        
        Parameters
        ----------
        k : int
            Number of taxa to select.
        
        Returns
        -------
        tuple[float, set[Any]]
            Tuple of (diversity_value, solution_set_of_node_ids).
        """
        self._fill_dp_table(k)
        
        root = self.network.root_node
        pd = self.table[root][k][frozenset()][0]
        
        solution = self._backtrack_solution(root, k, frozenset())
        
        return pd, solution


class AllPathsDiversity:
    """
    All-paths diversity measure.
    
    This measure computes the sum of branch lengths in the subnetwork
    induced by a set of taxa (all edges on paths from root to those taxa).
    """

    def compute_diversity(
        self,
        network: DirectedPhyNetwork,
        taxa: Set[str],
        **kwargs: Any,
    ) -> float:
        """
        Compute all-paths diversity of a set of taxa.
        
        The diversity is the sum of branch lengths of all edges in the
        subnetwork induced by the taxa (i.e., all edges on paths from
        root to the taxa).
        
        Parameters
        ----------
        network : DirectedPhyNetwork
            The phylogenetic network.
        taxa : Set[str]
            Set of taxa to compute diversity for.
        **kwargs : Any
            Additional parameters (currently unused).
        
        Returns
        -------
        float
            The all-paths diversity value (sum of branch lengths).
        
        Examples
        --------
        >>> from phylozoo.panda.measure import all_paths
        >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 3, 'v': 1, 'branch_length': 0.5},
        ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
        ...     ],
        ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        ... )
        >>> all_paths.compute_diversity(net, {"A", "B"})
        0.8
        """
        if not taxa:
            return 0.0
        
        # Get node IDs for the taxa
        # Note: Validation that taxa are in network is done in base.diversity()
        leaf_nodes: Set[Any] = set()
        for taxon in taxa:
            node_id = network.get_node_id(taxon)
            if node_id is None:
                raise ValueError(f"Taxon '{taxon}' not found in network")
            leaf_nodes.add(node_id)
        
        # Get all ancestors of the taxa (and the taxa themselves)
        dag = network._graph._graph
        nodes_to_keep: Set[Any] = set(leaf_nodes)
        for leaf in leaf_nodes:
            nodes_to_keep.update(nx.ancestors(dag, leaf))
        
        # Sum branch lengths of edges in the induced subnetwork
        total_diversity = 0.0
        
        for u, v, key, data in network._graph.edges(keys=True, data=True):
            if u in nodes_to_keep and v in nodes_to_keep:
                branch_length = data.get('branch_length')
                if branch_length is not None:
                    total_diversity += branch_length
                else:
                    # Default to 1.0 if no branch length specified
                    total_diversity += 1.0
        
        return total_diversity

    def solve_maximization(
        self,
        network: DirectedPhyNetwork,
        k: int,
        tree_extension: str | TreeExtension = "optimal_XP",
        **kwargs: Any,
    ) -> tuple[float, Set[str]]:
        """
        Solve MAPPD (Maximum A Posteriori Phylogenetic Diversity) using scanwidth-based DP.
        
        This implements the optimal algorithm for all-paths diversity using
        dynamic programming over a tree extension of the network.
        
        Parameters
        ----------
        network : DirectedPhyNetwork
            The phylogenetic network.
        k : int
            Number of taxa to select.
        tree_extension : str | TreeExtension, optional
            Method for computing tree extension or a precomputed TreeExtension.
            Currently only "optimal_XP" is supported (other methods not yet implemented).
            By default "optimal_XP".
        **kwargs : Any
            Additional parameters for tree extension computation.
        
        Returns
        -------
        tuple[float, Set[str]]
            Tuple of (diversity_value, solution_set_of_taxa).
        
        Raises
        ------
        PhyloZooValueError
            If k is invalid, network has parallel edges, network has 2-blobs,
            or network is invalid.
        PhyloZooRuntimeError
            If tree extension computation fails.
        PhyloZooNotImplementedError
            If tree_extension method is not "optimal_XP" (other methods not yet available).
        
        Notes
        -----
        Non-binary networks are automatically converted to binary using binary_resolution.
        If a precomputed TreeExtension is provided for a non-binary network, it will be
        recomputed after binary resolution.
        """
        if k < 0 or k > len(network.taxa):
            raise PhyloZooValueError(f"k must be between 0 and {len(network.taxa)}, got {k}")
        
        # Check for parallel edges
        if has_parallel_edges(network):
            raise PhyloZooValueError(
                "MAPPD algorithm cannot be applied to networks with parallel edges"
            )
        
        # Check for 2-blobs
        two_blobs = k_blobs(network, k=2, trivial=False, leaves=False)
        if two_blobs:
            raise PhyloZooValueError(
                f"MAPPD algorithm cannot be applied to networks with 2-blobs "
                f"(found {len(two_blobs)} 2-blob(s))"
            )
        
        # Apply binary resolution if network is not binary
        working_network = network
        if not is_binary(network):
            working_network = binary_resolution(network)
            # If network was non-binary, we can't use a precomputed tree extension
            if isinstance(tree_extension, TreeExtension):
                tree_extension = "optimal_XP"
        
        # Compute tree extension if needed
        if not isinstance(tree_extension, TreeExtension):
            if tree_extension != "optimal_XP":
                raise PhyloZooNotImplementedError(
                    f"Tree extension method '{tree_extension}' not yet implemented. "
                    "Only 'optimal_XP' is currently supported."
                )
            
            dag = ScanwidthDAG(working_network._graph._graph)
            res = dag.optimal_scanwidth(**kwargs)
            
            if res[0] is None:
                raise PhyloZooRuntimeError("Failed to compute tree extension")
            
            scanwidth, extension = res
            tree_extension = extension.canonical_tree_extension()
        
        # Create DP instance and solve
        dp_instance = _DPInstance(working_network, tree_extension)
        pd_value, solution_node_ids = dp_instance.solve(k)
        
        # Convert node IDs to taxa labels
        solution_taxa: Set[str] = set()
        for node_id in solution_node_ids:
            label = working_network.get_label(node_id)
            if label is not None:
                solution_taxa.add(label)
            else:
                # Should not happen for leaves, but handle gracefully
                raise PhyloZooRuntimeError(f"Leaf node {node_id} has no label")
        
        return pd_value, solution_taxa


# Convenience instance
all_paths = AllPathsDiversity()

