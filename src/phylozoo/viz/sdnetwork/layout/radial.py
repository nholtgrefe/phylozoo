"""
Radial layout algorithm for SemiDirectedPhyNetwork trees.

This module implements a circular/radial layout for phylogenetic trees.
Nodes are positioned in a circular arrangement with the root at the center
and leaves on the outer circle.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, TypeVar

import networkx as nx

from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.core.network.sdnetwork.derivations import to_d_network
from phylozoo.utils.exceptions import (
    PhyloZooLayoutError,
    PhyloZooValueError,
)

from .base import RadialLayout
from ..rendering.routes import compute_radial_routes

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

T = TypeVar('T')


def compute_radial_layout(
    network: 'SemiDirectedPhyNetwork',
    radius: float = 1.0,
    start_angle: float = 0.0,
    angle_direction: str = 'clockwise',
) -> RadialLayout:
    """
    Compute a radial (circular) layout for a SemiDirectedPhyNetwork tree.

    This function positions nodes in a circular arrangement:
    - Root at the center (radius 0)
    - Leaves on the outer circle
    - Internal nodes on concentric circles based on their depth from root

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to layout. Must be a tree.
    radius : float, optional
        Maximum radius for leaf nodes. By default 1.0.
    start_angle : float, optional
        Starting angle in radians for the first leaf. By default 0.0.
    angle_direction : str, optional
        Direction for angle progression: 'clockwise' or 'counterclockwise'.
        By default 'clockwise'.

    Returns
    -------
    RadialLayout
        The computed radial layout.

    Raises
    ------
    PhyloZooLayoutError
        If network is empty, not a tree, or layout computation fails.
    PhyloZooValueError
        If angle_direction is invalid.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork.layout import compute_radial_layout
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_radial_layout(net)
    >>> len(layout.positions)
    3
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")

    # Check if network is a tree
    if not is_tree(network):
        raise PhyloZooLayoutError(
            "Radial layout is only supported for trees. "
            "Use is_tree() to check if the network is a tree."
        )

    if angle_direction not in ('clockwise', 'counterclockwise'):
        raise PhyloZooValueError(
            f"angle_direction must be 'clockwise' or 'counterclockwise', got '{angle_direction}'"
        )

    # Convert to directed network to get root and tree structure
    # Suppress warnings from intermediate networks during conversion
    import warnings
    from phylozoo.utils.exceptions import (
        PhyloZooEmptyNetworkWarning,
        PhyloZooSingleNodeNetworkWarning,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=PhyloZooEmptyNetworkWarning)
        warnings.filterwarnings('ignore', category=PhyloZooSingleNodeNetworkWarning)
        d_network = to_d_network(network)
    root = d_network.root_node

    # Build tree structure from directed network
    G = nx.DiGraph()
    for node in d_network._graph.nodes:
        G.add_node(node)
    for u, v, key in d_network._graph.edges(keys=True):
        G.add_edge(u, v)

    # Compute depths from root
    depths: dict[T, int] = {}
    max_depth = 0

    def compute_depth(node: T, depth: int) -> None:
        """Recursively compute depth of nodes."""
        depths[node] = depth
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in G.successors(node):
            compute_depth(child, depth + 1)

    compute_depth(root, 0)

    # Get leaves in order (for angle assignment)
    leaves = list(d_network.leaves)
    num_leaves = len(leaves)

    if num_leaves == 0:
        # Single node network
        positions: dict[T, tuple[float, float]] = {root: (0.0, 0.0)}
        edge_routes: dict[tuple[T, T, int], 'EdgeRoute'] = {}
        return RadialLayout(
            network=network,
            positions=positions,
            edge_routes=edge_routes,
            algorithm='radial',
            parameters={
                'radius': radius,
                'start_angle': start_angle,
                'angle_direction': angle_direction,
            },
        )

    # Compute angles for leaves (evenly distributed around circle)
    angle_step = 2 * math.pi / num_leaves
    direction_mult = 1.0 if angle_direction == 'clockwise' else -1.0

    # Assign angles to leaves based on their subtree order
    # Use a recursive approach to assign angles based on subtree structure
    leaf_angles: dict[T, float] = {}
    leaf_index = 0

    def assign_leaf_angles(node: T) -> None:
        """Recursively assign angles to leaves based on subtree structure."""
        nonlocal leaf_index
        children = list(G.successors(node))
        if not children:
            # Leaf node
            angle = start_angle + leaf_index * angle_step * direction_mult
            leaf_angles[node] = angle
            leaf_index += 1
        else:
            # Internal node - process children in order
            for child in children:
                assign_leaf_angles(child)

    assign_leaf_angles(root)

    # Compute positions for all nodes
    positions: dict[T, tuple[float, float]] = {}

    def compute_positions(node: T) -> None:
        """Recursively compute positions for nodes."""
        if node in positions:
            return

        children = list(G.successors(node))
        if not children:
            # Leaf node
            angle = leaf_angles[node]
            node_radius = radius
            x = node_radius * math.cos(angle)
            y = node_radius * math.sin(angle)
            positions[node] = (x, y)
        else:
            # Internal node - first compute positions of all children
            for child in children:
                compute_positions(child)
            
            # Then position this node based on children
            child_positions = [positions[child] for child in children]
            
            # Compute average angle of children
            child_angles = [
                math.atan2(y, x) for x, y in child_positions
            ]
            # Normalize angles to [0, 2π)
            child_angles = [a if a >= 0 else a + 2 * math.pi for a in child_angles]
            avg_angle = sum(child_angles) / len(child_angles)
            
            # Position at radius proportional to depth
            # Root (depth 0) stays at center
            if node == root:
                positions[node] = (0.0, 0.0)
            elif max_depth > 0:
                node_radius = (depths[node] / max_depth) * radius
                x = node_radius * math.cos(avg_angle)
                y = node_radius * math.sin(avg_angle)
                positions[node] = (x, y)
            else:
                positions[node] = (0.0, 0.0)

    # Compute positions for all nodes starting from root
    compute_positions(root)
    
    # Filter positions to only include nodes from the original network
    # (to_d_network may create subdivision nodes that we don't want to position)
    original_nodes = set(network._graph.nodes)
    filtered_positions: dict[T, tuple[float, float]] = {
        node: pos for node, pos in positions.items() if node in original_nodes
    }

    # Compute edge routes (using filtered positions)
    edge_routes = compute_radial_routes(network, filtered_positions)

    return RadialLayout(
        network=network,
        positions=filtered_positions,
        edge_routes=edge_routes,
        algorithm='radial',
        parameters={
            'radius': radius,
            'start_angle': start_angle,
            'angle_direction': angle_direction,
        },
    )
