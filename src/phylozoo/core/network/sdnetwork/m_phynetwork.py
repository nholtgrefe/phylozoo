"""
Mixed network module.

This module provides classes and functions for working with mixed phylogenetic networks.
"""

from typing import List, Optional, Set, Tuple

from ..dnetwork.d_phynetwork import DirectedPhyNetwork


class MixedPhyNetwork(DirectedPhyNetwork):
    """
    A mixed phylogenetic network.
    
    This class is immutable after initialization. To create or modify a network,
    build it using MixedMultiGraph and then create a MixedPhyNetwork from it.
    """

    def __init__(
        self,
        edges: Optional[List[Tuple[str, str]]] = None,
        taxa: Optional[dict] = None,
    ) -> None:
        """
        Initialize a mixed network.

        Parameters
        ----------
        edges : Optional[List[Tuple[str, str]]], optional
            List of edges as tuples, by default None
        taxa : Optional[dict], optional
            Taxon labels for leaves, by default None
        """
        # Convert None to empty list for parent class requirement
        edges_list = edges if edges is not None else []
        super().__init__(edges=edges_list, taxa=taxa)
        # TODO: Add mixed network specific attributes

    def __repr__(self) -> str:
        """
        Return string representation of the network.

        Returns
        -------
        str
            String representation
        """
        return (
            f"MixedPhyNetwork(nodes={self.number_of_nodes()}, "
            f"edges={self.number_of_edges()})"
        )

