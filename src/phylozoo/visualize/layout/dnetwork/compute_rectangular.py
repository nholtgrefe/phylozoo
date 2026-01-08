"""
Computation functions for rectangular layouts of DirectedPhyNetwork.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

import networkx as nx

from ..base import EdgeRoute, EdgeType
from .base import RectangularDNetLayout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


def compute_rectangular_dnet_layout(
    network: 'DirectedPhyNetwork',
    layer_gap: float = 1.5,
    node_gap: float = 1.0,
    iterations: int = 50,
    seed: int | None = None,
    orientation: str = 'top-bottom',
    curve_strength: float = 0.3,
    node_width: float = 0.1,
) -> RectangularDNetLayout:
    """
    Compute rectangular layout for a DirectedPhyNetwork.

    This function computes node positions and edge routing using the
    rectangular layout algorithm. It positions nodes in strict horizontal
    layers with crossing minimization.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    layer_gap : float, optional
        Vertical spacing between layers. By default 1.5.
    node_gap : float, optional
        Horizontal spacing between nodes within a layer. By default 1.0.
    iterations : int, optional
        Number of iterations for crossing minimization. By default 50.
    seed : int | None, optional
        Random seed for reproducibility. By default None.
    orientation : str, optional
        Layout direction: 'top-bottom' or 'left-right'. By default 'top-bottom'.
    curve_strength : float, optional
        Strength of curve for hybrid edges (0.0 = straight, 1.0 = strong).
        By default 0.3.
    node_width : float, optional
        Width of nodes for computing hybrid edge offsets. By default 0.1.

    Returns
    -------
    RectangularDNetLayout
        The computed layout.

    Raises
    ------
    ValueError
        If network is not a DAG or has no root node.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize.layout.dnetwork import compute_rectangular_dnet_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_rectangular_dnet_layout(net)
    >>> len(layout.positions)
    3
    """
    # Step 1: Compute node positions (reuse existing function)
    # Note: rectangular_layout moved to sandbox, import from there if needed
    import sys
    from pathlib import Path
    # Get project root (go up from src/phylozoo/visualize/layout/dnetwork)
    project_root = Path(__file__).parent.parent.parent.parent.parent.parent
    sandbox_path = project_root / 'sandbox'
    if str(sandbox_path) not in sys.path:
        sys.path.insert(0, str(sandbox_path))
    from rectangular_layout import compute_rectangular_layout

    positions = compute_rectangular_layout(
        network,
        layer_gap=layer_gap,
        node_gap=node_gap,
        iterations=iterations,
        seed=seed,
        orientation=orientation,
        use_branch_lengths=False,  # Not yet implemented
    )

    # Step 2: Compute layers (needed for edge routing)
    G_nx = nx.DiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        G_nx.add_edge(u, v)

    root = network.root_node
    layers = _assign_layers(G_nx, root)

    # Step 3: Compute edge routing
    edge_routes = _compute_edge_routes(
        network, positions, layers, curve_strength, node_width, orientation
    )

    # Step 4: Create layout object
    return RectangularDNetLayout(
        graph=network,  # Pass as graph to parent
        positions=positions,
        edge_routes=edge_routes,
        layers=layers,
        orientation=orientation,
        algorithm='rectangular',
        parameters={
            'layer_gap': layer_gap,
            'node_gap': node_gap,
            'iterations': iterations,
            'seed': seed,
            'curve_strength': curve_strength,
            'node_width': node_width,
        },
    )


# ============================================================================
# Helper Functions
# ============================================================================


def _assign_layers(G_nx: nx.DiGraph, root: T) -> dict[T, int]:
    """
    Assign nodes to layers based on longest path from root.

    Parameters
    ----------
    G_nx : nx.DiGraph
        NetworkX directed graph representation.
    root : T
        Root node of the network.

    Returns
    -------
    dict[T, int]
        Dictionary mapping node to layer number (0 = root layer).
    """
    layers: dict[T, int] = {root: 0}

    # Use topological sort to process nodes in correct order
    topo_order = list(nx.topological_sort(G_nx))

    # Assign layers based on longest path
    for node in topo_order:
        if node == root:
            continue

        # Get all predecessors
        predecessors = list(G_nx.predecessors(node))
        if not predecessors:
            layers[node] = 0
        else:
            # Layer is max of all predecessor layers + 1
            max_pred_layer = max(layers.get(p, 0) for p in predecessors)
            layers[node] = max_pred_layer + 1

    # Move all leaves to the bottom layer
    leaves = [n for n in G_nx.nodes if G_nx.out_degree(n) == 0]
    max_layer = max(layers.values()) if layers else 0

    for leaf in leaves:
        layers[leaf] = max_layer

    # Recompute layers for nodes that are now "above" leaves
    changed = True
    while changed:
        changed = False
        for node in topo_order:
            if node == root or node in leaves:
                continue

            predecessors = list(G_nx.predecessors(node))
            if predecessors:
                max_pred_layer = max(layers.get(p, 0) for p in predecessors)
                new_layer = max_pred_layer + 1
                if new_layer >= max_layer:
                    new_layer = max_layer - 1
                if layers.get(node, 0) != new_layer:
                    layers[node] = new_layer
                    changed = True

    return layers


def _compute_edge_routes(
    network: 'DirectedPhyNetwork',
    positions: dict[T, tuple[float, float]],
    layers: dict[T, int],
    curve_strength: float,
    node_width: float,
    orientation: str,
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Compute edge routing for all edges in the network.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    positions : dict[T, tuple[float, float]]
        Node positions.
    layers : dict[T, int]
        Layer assignments.
    curve_strength : float
        Strength of curve for hybrid edges.
    node_width : float
        Width of nodes for computing offsets.
    orientation : str
        Layout orientation.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information.
    """
    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}
    hybrid_nodes = network.hybrid_nodes

    # First pass: count parallel edges to determine is_parallel
    parallel_counts: dict[tuple[T, T], int] = {}
    for u, v, key in network._graph.edges(keys=True):
        edge_pair = (u, v)
        parallel_counts[edge_pair] = parallel_counts.get(edge_pair, 0) + 1

    # Second pass: compute routes with edge type information
    for u, v, key in network._graph.edges(keys=True):
        is_hybrid = v in hybrid_nodes
        is_parallel = parallel_counts[(u, v)] > 1

        if is_hybrid:
            # Curved inlet for hybrid edge
            route = _compute_curved_inlet(
                positions[u],
                positions[v],
                layers.get(u, 0),
                layers.get(v, 0),
                curve_strength,
                node_width,
                orientation,
            )
            edge_routes[(u, v, key)] = EdgeRoute(
                edge_type=EdgeType(
                    is_directed=True,
                    is_hybrid=True,
                    is_parallel=is_parallel,
                ),
                points=route['points'],
                curve_control=route['control'],
            )
        else:
            # Rectangular routing for tree edge
            route_points = _compute_rectangular_route(
                positions[u],
                positions[v],
                layers.get(u, 0),
                layers.get(v, 0),
                orientation,
            )
            edge_routes[(u, v, key)] = EdgeRoute(
                edge_type=EdgeType(
                    is_directed=True,
                    is_hybrid=False,
                    is_parallel=is_parallel,
                ),
                points=route_points,
                curve_control=None,
            )

    return edge_routes


def _compute_rectangular_route(
    parent_pos: tuple[float, float],
    child_pos: tuple[float, float],
    parent_layer: int,
    child_layer: int,
    orientation: str,
) -> tuple[tuple[float, float], ...]:
    """
    Compute rectangular (orthogonal) route for a tree edge.

    Parameters
    ----------
    parent_pos : tuple[float, float]
        Parent node position.
    child_pos : tuple[float, float]
        Child node position.
    parent_layer : int
        Parent layer number.
    child_layer : int
        Child layer number.
    orientation : str
        Layout orientation.

    Returns
    -------
    tuple[tuple[float, float], ...]
        Polyline points for the route.
    """
    px, py = parent_pos
    cx, cy = child_pos

    if orientation == 'top-bottom':
        # Vertical then horizontal then vertical
        mid_y = (py + cy) / 2
        return ((px, py), (px, mid_y), (cx, mid_y), (cx, cy))
    else:  # left-right
        # Horizontal then vertical then horizontal
        mid_x = (px + cx) / 2
        return ((px, py), (mid_x, py), (mid_x, cy), (cx, cy))


def _compute_curved_inlet(
    parent_pos: tuple[float, float],
    hybrid_pos: tuple[float, float],
    parent_layer: int,
    hybrid_layer: int,
    curve_strength: float,
    node_width: float,
    orientation: str,
) -> dict[str, Any]:
    """
    Compute curved inlet route for a hybrid edge.

    Parameters
    ----------
    parent_pos : tuple[float, float]
        Parent node position.
    hybrid_pos : tuple[float, float]
        Hybrid node position.
    parent_layer : int
        Parent layer number.
    hybrid_layer : int
        Hybrid layer number.
    curve_strength : float
        Strength of the curve (0.0 = straight, 1.0 = strong).
    node_width : float
        Width of nodes for computing offsets.
    orientation : str
        Layout orientation.

    Returns
    -------
    dict[str, Any]
        Dictionary with 'points' and 'control' keys.
    """
    px, py = parent_pos
    hx, hy = hybrid_pos

    if orientation == 'top-bottom':
        # Compute offset for multiple parents (simple version)
        # In a full implementation, we'd track all parents and offset accordingly
        offset_x = curve_strength * node_width

        # Create curved path: start from parent, curve toward hybrid
        mid_y = (py + hy) / 2
        control_x = (px + hx) / 2 + offset_x
        control_y = mid_y

        points = (
            (px, py),
            (px, mid_y),
            (control_x, control_y),
            (hx, mid_y),
            (hx, hy),
        )

        return {
            'points': points,
            'control': (control_x, control_y),
        }
    else:  # left-right
        # Similar logic but swapped coordinates
        offset_y = curve_strength * node_width
        mid_x = (px + hx) / 2
        control_x = mid_x
        control_y = (py + hy) / 2 + offset_y

        points = (
            (px, py),
            (mid_x, py),
            (control_x, control_y),
            (mid_x, hy),
            (hx, hy),
        )

        return {
            'points': points,
            'control': (control_x, control_y),
        }

