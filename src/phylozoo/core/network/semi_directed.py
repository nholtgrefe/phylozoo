"""
Semi-directed network module.

This module provides classes and functions for working with semi-directed phylogenetic networks.
"""

from typing import List, Optional, Set, Tuple

from .directed import DirectedNetwork, validates


class SemiDirectedNetwork(DirectedNetwork):
    """
    A semi-directed phylogenetic network.
    
    Inherits validation context manager from DirectedNetwork.
    """

    def __init__(
        self, nodes: Optional[Set[str]] = None, edges: Optional[List[Tuple[str, str]]] = None
    ) -> None:
        """
        Initialize a semi-directed network.

        Parameters
        ----------
        nodes : Optional[Set[str]], optional
            Set of node identifiers, by default None
        edges : Optional[List[Tuple[str, str]]], optional
            List of edges as tuples, by default None
        """
        super().__init__(nodes, edges)
        self.undirected_edges: List[Tuple[str, str]] = []

    @validates
    def add_undirected_edge(self, node1: str, node2: str) -> None:
        """
        Add an undirected edge to the network.

        Parameters
        ----------
        node1 : str
            First node identifier
        node2 : str
            Second node identifier
        """
        if node1 not in self.nodes:
            self.add_node(node1)
        if node2 not in self.nodes:
            self.add_node(node2)
        self.undirected_edges.append((node1, node2))


    def _validate_network(self) -> None:
        """
        Validate entire semi-directed network state.
        
        Checks that the network is a valid semi-directed network (e.g., all edges
        reference existing nodes, network invariants are maintained, etc.).
        
        Raises
        ------
        ValueError
            If network is invalid
        """
        # Validate directed edges (parent validation)
        super()._validate_network()
        
        # Validate undirected edges
        errors = []
        for node1, node2 in self.undirected_edges:
            if node1 not in self.nodes:
                errors.append(f"Undirected edge node '{node1}' not in nodes")
            if node2 not in self.nodes:
                errors.append(f"Undirected edge node '{node2}' not in nodes")
        
        # Add other semi-directed network-specific validations as needed
        # - Check for self-loops in undirected edges if not allowed
        # - Check network topology invariants
        # - etc.
        
        if errors:
            raise ValueError(f"Invalid semi-directed network:\n" + "\n".join(errors))

    def __repr__(self) -> str:
        """
        Return string representation of the network.

        Returns
        -------
        str
            String representation
        """
        return (
            f"SemiDirectedNetwork(nodes={len(self.nodes)}, "
            f"directed_edges={len(self.edges)}, "
            f"undirected_edges={len(self.undirected_edges)})"
        )


def random_semi_directed_network(n_leaves: int, seed: Optional[int] = None) -> SemiDirectedNetwork:
    """
    Generate a random semi-directed network.

    Parameters
    ----------
    n_leaves : int
        Number of leaves in the network
    seed : Optional[int], optional
        Random seed, by default None

    Returns
    -------
    SemiDirectedNetwork
        A random semi-directed network (placeholder - returns empty network)

    Notes
    -----
    This is a placeholder function. Implement actual random generation logic here.
    """
    # Placeholder implementation
    return SemiDirectedNetwork()
