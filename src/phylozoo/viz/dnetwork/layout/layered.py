"""
Layered (Sugiyama-style) layout for DirectedPhyNetwork (pz-layered).

A self-contained equivalent of Graphviz ``dot``: nodes are assigned to layers
by their longest path from the root, edges spanning several layers are routed
through dummy nodes so they bend around intermediate layers, the order within
each layer is chosen by barycenter sweeps with adjacent swaps to reduce edge
crossings, and coordinates are assigned so that edges are as straight as
possible. Unlike ``pz-cladogram`` there is no tree backbone: nodes of a layer are
ordered freely, which gives fewer crossings on heavily reticulate networks at
the cost of sibling subtrees not staying contiguous.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import networkx as nx

from phylozoo.utils.exceptions import (
    PhyloZooLayoutError,
    PhyloZooNetworkStructureError,
    PhyloZooValueError,
)
from phylozoo.viz._layout_utils import normalize_positions
from phylozoo.viz._types import EdgeRoute, EdgeType

from .base import DNetLayout

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

Layers = list[list[Any]]


class _Dummy:
    """Bend point of an edge that spans several layers."""

    __slots__ = ()


def compute_pz_layered_layout(
    network: "DirectedPhyNetwork",
    layer_gap: float = 1.0,
    node_gap: float = 1.0,
    sweeps: int = 8,
    direction: str = "TD",
    align_leaves: bool = False,
    x_scale: float = 1.0,
    y_scale: float = 1.0,
    rectangular: bool = False,
) -> DNetLayout:
    """
    Layout a phylogenetic network (DAG) with the layered (Sugiyama) method of ``dot``.

    This is a custom PhyloZoo layout algorithm (pz-layered).

    The algorithm:

    1. Assigns every node to a layer, its longest path from the root, so every
       edge points to a strictly lower layer and hybrid nodes sit right below
       their lowest parent.
    2. Replaces every edge that spans more than one layer by a chain through
       dummy nodes, one per intermediate layer, so long edges are ordered and
       routed like nodes instead of cutting across.
    3. Orders each layer by alternating downward and upward barycenter sweeps
       (a node moves to the mean position of its neighbours in the layer just
       processed), each followed by swapping adjacent nodes while that lowers
       the number of crossings with the neighbouring layers. The ordering with
       the fewest crossings over all sweeps is kept.
    4. Assigns coordinates within layers: a few passes pull every node towards
       the mean of its neighbours (dummy nodes towards their edge, keeping it
       straight) while keeping nodes ``node_gap`` apart.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    layer_gap : float, optional
        Spacing between layers. By default 1.0.
    node_gap : float, optional
        Minimum spacing between nodes of the same layer. By default 1.0.
    sweeps : int, optional
        Number of down-and-up ordering sweeps. By default 8.
    direction : str, optional
        Layout direction: 'TD' (root at the top) or 'LR' (root on the left).
        By default 'TD'.
    align_leaves : bool, optional
        Put all leaves on the bottom layer. By default False (leaves stay at
        their own depth, as in ``dot``).
    x_scale : float, optional
        Scaling factor applied to x coordinates. By default 1.0.
    y_scale : float, optional
        Scaling factor applied to y coordinates. By default 1.0.
    rectangular : bool, optional
        Draw every edge segment as an orthogonal elbow: along the layer of its
        upper end and then down; the last segment into a hybrid node instead
        runs down first and enters the hybrid sideways. By default False.

    Returns
    -------
    DNetLayout
        The computed layout, with polyline edge routes and
        ``algorithm='pz-layered'``. Positions are centred and scaled to unit range.

    Raises
    ------
    PhyloZooLayoutError
        If network is empty.
    PhyloZooNetworkStructureError
        If network is not a DAG or has no root node.
    PhyloZooValueError
        If direction is invalid.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz.dnetwork.layout import compute_pz_layered_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(5, 3), (5, 4), (3, 2), (4, 2), (2, 1), (3, 6), (4, 7)],
    ...     nodes=[(1, {'label': 'A'}), (6, {'label': 'B'}), (7, {'label': 'C'})]
    ... )
    >>> layout = compute_pz_layered_layout(net)
    >>> layout.algorithm
    'pz-layered'
    >>> len(layout.positions) == net.number_of_nodes()
    True
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")
    if direction.upper() not in ("TD", "LR"):
        raise PhyloZooValueError(f"direction must be 'TD' or 'LR', got '{direction}'")

    graph: nx.DiGraph = nx.DiGraph()
    graph.add_nodes_from(network._graph.nodes)
    graph.add_edges_from((u, v) for u, v in network._graph.edges if u != v)
    if not nx.is_directed_acyclic_graph(graph):
        raise PhyloZooNetworkStructureError("Network must be a DAG")
    roots = [n for n in graph if graph.in_degree(n) == 0]
    if not roots:
        raise PhyloZooNetworkStructureError("No root found (no node with in-degree 0).")

    # --- 1. Layers ---
    topo = list(nx.topological_sort(graph))
    depth: dict[Any, int] = {}
    for n in topo:
        depth[n] = max((depth[p] + 1 for p in graph.predecessors(n)), default=0)
    if align_leaves:
        bottom = max(depth.values())
        for n in graph:
            if graph.out_degree(n) == 0:
                depth[n] = bottom

    # --- 2. Proper layered graph with dummy nodes ---
    proper: nx.DiGraph = nx.DiGraph()
    proper.add_nodes_from(graph)
    chain: dict[tuple[Any, Any], list[Any]] = {}
    for u, v in graph.edges:
        bends = [_Dummy() for _ in range(depth[v] - depth[u] - 1)]
        for i, d in enumerate(bends):
            depth[d] = depth[u] + i + 1
        chain[(u, v)] = bends
        nx.add_path(proper, [u, *bends, v])

    layers: Layers = [[] for _ in range(max(depth.values()) + 1)]
    for n in nx.dfs_preorder_nodes(proper, roots[0]):
        layers[depth[n]].append(n)
    for n in proper:  # nodes unreachable from the first root (extra roots)
        if n not in layers[depth[n]]:
            layers[depth[n]].append(n)

    # --- 3. Ordering ---
    best_layers = [list(layer) for layer in layers]
    best_crossings = _total_crossings(proper, layers)
    for _ in range(sweeps):
        for downward in (True, False):
            _barycenter_sweep(proper, layers, downward)
            _transpose(proper, layers)
            crossings = _total_crossings(proper, layers)
            if crossings < best_crossings:
                best_crossings = crossings
                best_layers = [list(layer) for layer in layers]
    layers = best_layers

    # --- 4. Coordinates ---
    x = _assign_coordinates(proper, layers, node_gap)
    max_depth = len(layers) - 1
    all_positions: dict[Any, tuple[float, float]] = {}
    for n in proper:
        if direction.upper() == "TD":
            all_positions[n] = (x[n] * x_scale, (max_depth - depth[n]) * layer_gap * y_scale)
        else:
            all_positions[n] = (depth[n] * layer_gap * x_scale, -x[n] * y_scale)
    all_positions = normalize_positions(all_positions)
    positions = {n: all_positions[n] for n in graph}

    # --- Routes through the dummy nodes ---
    counts: dict[tuple[Any, Any], int] = {}
    for u, v, _ in network._graph.edges(keys=True):
        counts[(u, v)] = counts.get((u, v), 0) + 1
    routes: dict[tuple[Any, Any, int], EdgeRoute] = {}
    backbone_edges: set[tuple[Any, Any, int]] = set()
    reticulate_edges: set[tuple[Any, Any, int]] = set()
    for u, v, key in network._graph.edges(keys=True):
        is_hybrid = v in network.hybrid_nodes
        (reticulate_edges if is_hybrid else backbone_edges).add((u, v, key))
        points = tuple(all_positions[n] for n in (u, *chain.get((u, v), []), v))
        if rectangular:
            points = _orthogonal(points, is_hybrid, direction.upper() == "TD")
        routes[(u, v, key)] = EdgeRoute(
            edge_type=EdgeType(
                is_directed=True, is_hybrid=is_hybrid, is_parallel=counts[(u, v)] > 1
            ),
            points=points,
        )

    return DNetLayout(
        network=network,
        positions=positions,
        edge_routes=routes,
        backbone_edges=backbone_edges,
        reticulate_edges=reticulate_edges,
        algorithm="pz-layered",
        parameters={
            "layer_gap": layer_gap,
            "node_gap": node_gap,
            "sweeps": sweeps,
            "direction": direction,
            "align_leaves": align_leaves,
            "x_scale": x_scale,
            "y_scale": y_scale,
            "rectangular": rectangular,
        },
    )


def _orthogonal(
    points: tuple[tuple[float, float], ...], into_hybrid: bool, top_down: bool
) -> tuple[tuple[float, float], ...]:
    """Turn a polyline into elbows; the last segment into a hybrid arrives sideways."""
    route: list[tuple[float, float]] = [points[0]]
    for k, (a, b) in enumerate(zip(points, points[1:])):
        last = k == len(points) - 2
        if into_hybrid and last:
            corner = (a[0], b[1]) if top_down else (b[0], a[1])
        else:
            corner = (b[0], a[1]) if top_down else (a[0], b[1])
        if corner != a and corner != b:
            route.append(corner)
        route.append(b)
    return tuple(route)


def _crossings_between(proper: nx.DiGraph, upper: list[Any], lower: list[Any]) -> int:
    """Number of crossings among the edges between two consecutive layers."""
    pos_up = {n: i for i, n in enumerate(upper)}
    pos_low = {n: i for i, n in enumerate(lower)}
    pairs = sorted(
        (pos_up[u], pos_low[v]) for u in upper for v in proper.successors(u) if v in pos_low
    )
    crossings = 0
    for i, (_, b) in enumerate(pairs):
        for _, d in pairs[i + 1 :]:
            if d < b:
                crossings += 1
    return crossings


def _total_crossings(proper: nx.DiGraph, layers: Layers) -> int:
    return sum(_crossings_between(proper, layers[i], layers[i + 1]) for i in range(len(layers) - 1))


def _barycenter_sweep(proper: nx.DiGraph, layers: Layers, downward: bool) -> None:
    """Reorder every layer by the mean position of its neighbours in the previous layer."""
    indices = range(1, len(layers)) if downward else range(len(layers) - 2, -1, -1)
    for i in indices:
        reference = {n: k for k, n in enumerate(layers[i - 1 if downward else i + 1])}
        neighbours = proper.predecessors if downward else proper.successors
        current = {n: k for k, n in enumerate(layers[i])}

        def barycenter(n: Any) -> float:
            ref = [reference[m] for m in neighbours(n) if m in reference]
            return sum(ref) / len(ref) if ref else float(current[n])

        layers[i].sort(key=lambda n: (barycenter(n), current[n]))


def _transpose(proper: nx.DiGraph, layers: Layers) -> None:
    """Swap adjacent nodes within layers while that reduces crossings."""
    improved = True
    while improved:
        improved = False
        for i, layer in enumerate(layers):
            for k in range(len(layer) - 1):
                before = _local_crossings(proper, layers, i)
                layer[k], layer[k + 1] = layer[k + 1], layer[k]
                if _local_crossings(proper, layers, i) < before:
                    improved = True
                else:
                    layer[k], layer[k + 1] = layer[k + 1], layer[k]


def _local_crossings(proper: nx.DiGraph, layers: Layers, i: int) -> int:
    total = 0
    if i > 0:
        total += _crossings_between(proper, layers[i - 1], layers[i])
    if i + 1 < len(layers):
        total += _crossings_between(proper, layers[i], layers[i + 1])
    return total


def _assign_coordinates(proper: nx.DiGraph, layers: Layers, gap: float) -> dict[Any, float]:
    """Pull nodes towards their neighbours' mean while keeping them ``gap`` apart (priority method)."""
    x: dict[Any, float] = {}
    for layer in layers:
        for k, n in enumerate(layer):
            x[n] = k * gap

    def place(layer: list[Any], neighbours: Any) -> None:
        desired = {}
        for n in layer:
            if isinstance(n, _Dummy):  # keep the edge straight: mean of both chain neighbours
                ref = [x[m] for m in (*proper.predecessors(n), *proper.successors(n))]
            else:
                ref = [x[m] for m in neighbours(n) if m in x]
            desired[n] = sum(ref) / len(ref) if ref else x[n]
        left: dict[Any, float] = {}
        for k, n in enumerate(layer):
            left[n] = desired[n] if k == 0 else max(desired[n], left[layer[k - 1]] + gap)
        right: dict[Any, float] = {}
        for k in range(len(layer) - 1, -1, -1):
            n = layer[k]
            right[n] = (
                desired[n] if k == len(layer) - 1 else min(desired[n], right[layer[k + 1]] - gap)
            )
        for n in layer:
            x[n] = (left[n] + right[n]) / 2

    for _ in range(4):
        for layer in layers[1:]:
            place(layer, proper.predecessors)
        for layer in reversed(layers[:-1]):
            place(layer, proper.successors)
    return x
