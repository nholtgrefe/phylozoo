"""
DAG layout algorithm for DirectedPhyNetwork (pz-cladogram).

A layered drawing with a tree backbone: every node except the root hangs below
one of its parents (its lowest parent), and the remaining edges are drawn as
straight reticulate edges. The left-to-right order of the leaves is chosen by
a hybrid-aware barycenter heuristic followed by a local search on the exact
number of edge crossings, which keeps reticulate edges short and drawings
deterministic.
"""

from __future__ import annotations

import random
from typing import TYPE_CHECKING, Any, Callable, TypeVar

import networkx as nx
import numpy as np

from phylozoo.utils.exceptions import (
    PhyloZooLayoutError,
    PhyloZooNetworkStructureError,
    PhyloZooValueError,
)
from phylozoo.viz._layout_utils import count_crossings as _count_crossings
from phylozoo.viz._layout_utils import normalize_positions, sort_key

from .base import DNetLayout
from .routes import compute_backbone_routes, compute_hybrid_routes, compute_rectangular_routes

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar("T")
Children = dict[Any, list[Any]]

# Above this many edges the exact crossing-count local search is skipped.
_LOCAL_SEARCH_MAX_EDGES = 400


def compute_pz_cladogram_layout(
    network: "DirectedPhyNetwork",
    layer_gap: float = 1.0,
    leaf_gap: float = 1.0,
    trials: int = 10,
    seed: int | None = 0,
    direction: str = "TD",
    x_scale: float = 1.0,
    y_scale: float = 1.0,
    align_leaves: bool = True,
    rectangular: bool = True,
    horizontal_reticulations: bool | None = None,
) -> DNetLayout:
    """
    Layout a phylogenetic network (DAG) as a layered drawing with a tree backbone.

    This is a custom PhyloZoo layout algorithm (pz-cladogram).

    The algorithm:

    1. Assigns every node to a layer (its longest path from the root); leaves
       are moved to the bottom layer unless ``align_leaves`` is False.
    2. Hangs every node below its lowest parent, which gives a spanning tree
       (the backbone). Edges from the other parents of hybrid nodes are the
       reticulate edges.
    3. Orders the children of every backbone node: a few barycenter passes
       pull subtrees towards the partners of their reticulate edges, then a
       local search swaps adjacent siblings whenever that lowers the exact
       number of edge crossings (skipped for very large networks). ``trials``
       restarts of this procedure, each from a random order of the nodes that
       can affect crossings (ancestors of reticulate-edge endpoints), are
       compared and the best kept.
    4. Places leaves left to right, internal nodes above the mean of their
       children, and maps layers to the vertical (``'TD'``) or horizontal
       (``'LR'``) axis.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    layer_gap : float, optional
        Spacing between layers. By default 1.0.
    leaf_gap : float, optional
        Spacing between consecutive leaves. By default 1.0.
    trials : int, optional
        Number of ordering attempts; the first starts from the network's node
        order, the others from random orders of the crossing-relevant nodes.
        Trees need one. By default 10.
    seed : int | None, optional
        Random seed for the restarts; None draws a fresh seed. By default 0.
    direction : str, optional
        Layout direction: 'TD' (root at the top) or 'LR' (root on the left).
        By default 'TD'.
    x_scale : float, optional
        Scaling factor applied to x coordinates. By default 1.0.
    y_scale : float, optional
        Scaling factor applied to y coordinates. By default 1.0.
    align_leaves : bool, optional
        Put all leaves on the bottom layer. By default True.
    rectangular : bool, optional
        Draw edges as orthogonal elbows (backbone edges run along the parent's
        layer and then into the child, reticulate edges enter the hybrid node
        sideways along its layer); ``False`` draws straight lines. By default True.
    horizontal_reticulations : bool | None, optional
        Move the parent of every reticulate edge down to its hybrid's layer
        whenever the network allows it (its other descendants are pushed down
        as needed), so the reticulate edge is horizontal. None (default) means
        True for ``rectangular`` drawings and False otherwise.
    Returns
    -------
    DNetLayout
        The computed layout. Positions are centred at the origin and scaled to
        unit range.

    Raises
    ------
    PhyloZooLayoutError
        If network is empty or layout computation fails.
    PhyloZooNetworkStructureError
        If network is not a DAG or has no root node.
    PhyloZooValueError
        If direction is invalid (must be 'TD' or 'LR').

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz.dnetwork.layout import compute_pz_cladogram_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_pz_cladogram_layout(net)
    >>> len(layout.positions)
    3
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")
    if direction.upper() not in ("TD", "LR"):
        raise PhyloZooValueError(f"direction must be 'TD' or 'LR', got '{direction}'")

    graph: nx.DiGraph = nx.DiGraph()
    graph.add_nodes_from(sorted(network._graph.nodes, key=sort_key))
    graph.add_edges_from(
        sorted(
            ((u, v) for u, v in network._graph.edges if u != v),
            key=lambda e: tuple(map(sort_key, e)),
        )
    )
    if not nx.is_directed_acyclic_graph(graph):
        raise PhyloZooNetworkStructureError("Network must be a DAG")
    roots = [n for n in graph if graph.in_degree(n) == 0]
    if not roots:
        raise PhyloZooNetworkStructureError("No root found (no node with in-degree 0).")
    root = roots[0]

    # --- Layers and backbone ---
    topo = list(nx.topological_sort(graph))
    rank = {n: i for i, n in enumerate(topo)}
    depth: dict[Any, int] = {}
    for n in topo:
        depth[n] = max((depth[p] + 1 for p in graph.predecessors(n)), default=0)
    kids: Children = {n: [] for n in graph}
    for n in topo:
        parents = list(graph.predecessors(n))
        if parents:
            kids[max(parents, key=lambda p: (depth[p], -rank[p]))].append(n)
    for n in kids:
        kids[n].sort(key=rank.__getitem__)
    backbone_pairs = {(p, c) for p, cs in kids.items() for c in cs}
    edges = list(graph.edges)
    reticulate = [e for e in edges if e not in backbone_pairs]
    if horizontal_reticulations is None:
        horizontal_reticulations = rectangular
    if horizontal_reticulations:
        _level_reticulations(graph, depth, reticulate, topo)
    if align_leaves:
        bottom = max(depth.values())
        for n in graph:
            if graph.out_degree(n) == 0:
                depth[n] = bottom

    # --- Child ordering ---
    # Only nodes above an endpoint of a reticulate edge can change the number of
    # crossings; restarts shuffle those alone (a tree needs no restarts at all).
    parent_of = {c: p for p, cs in kids.items() for c in cs}
    relevant: set[Any] = set()
    for u, v in reticulate:
        for n in (u, v):
            while n in parent_of:
                n = parent_of[n]
                relevant.add(n)
    shuffle_nodes = [n for n in topo if n in relevant and len(kids[n]) > 1]
    rng = random.Random(seed)
    best_score: tuple[int, float] | None = None
    best_kids: Children = kids
    for trial in range(max(1, trials) if shuffle_nodes else 1):
        if trial > 0:
            for n in shuffle_nodes:
                rng.shuffle(kids[n])
        for _ in range(3):
            _barycenter_pass(root, kids, reticulate, leaf_gap)
        if len(edges) <= _LOCAL_SEARCH_MAX_EDGES:
            score = _local_search(
                root, kids, edges, reticulate, depth, leaf_gap, nodes=shuffle_nodes
            )
        else:
            score = _score(root, kids, edges, reticulate, depth, leaf_gap)
        if best_score is None or score < best_score:
            best_score = score
            best_kids = {n: list(cs) for n, cs in kids.items()}
    x = _x_coordinates(root, best_kids, leaf_gap)

    # --- Final coordinates ---
    max_depth = max(depth.values())
    positions: dict[Any, tuple[float, float]] = {}
    for n in graph:
        if direction.upper() == "TD":
            positions[n] = (x[n] * x_scale, (max_depth - depth[n]) * layer_gap * y_scale)
        else:
            positions[n] = (depth[n] * layer_gap * x_scale, -x[n] * y_scale)
    positions = normalize_positions(positions)

    # --- Edge routes ---
    backbone_edges: set[tuple[Any, Any, int]] = set()
    reticulate_edges: set[tuple[Any, Any, int]] = set()
    for u, v, key in network._graph.edges(keys=True):
        (backbone_edges if (u, v) in backbone_pairs else reticulate_edges).add((u, v, key))
    if rectangular:
        routes = compute_rectangular_routes(
            network, positions, backbone_edges, reticulate_edges, direction
        )
    else:
        routes = compute_backbone_routes(network, positions, backbone_edges)
        routes.update(compute_hybrid_routes(network, positions, reticulate_edges))

    return DNetLayout(
        network=network,
        positions=positions,
        edge_routes=routes,
        backbone_edges=backbone_edges,
        reticulate_edges=reticulate_edges,
        algorithm="pz-cladogram",
        parameters={
            "layer_gap": layer_gap,
            "leaf_gap": leaf_gap,
            "trials": trials,
            "seed": seed,
            "direction": direction,
            "x_scale": x_scale,
            "y_scale": y_scale,
            "align_leaves": align_leaves,
            "rectangular": rectangular,
            "horizontal_reticulations": horizontal_reticulations,
        },
    )


def _level_reticulations(
    graph: nx.DiGraph, depth: dict[Any, int], reticulate: list[tuple[Any, Any]], topo: list[Any]
) -> set[tuple[Any, Any]]:
    """
    Lower reticulate parents onto their hybrid's layer, in place, where the DAG allows it.

    Each reticulate edge is tried in turn (hybrids nearest the root first): the
    layering is recomputed with that edge constrained to zero layer difference
    and every other edge to at least one. A constraint set that cannot be
    satisfied (the parent reaches the hybrid by another path) is skipped.
    Returns the edges that ended up horizontal.
    """
    horizontal: set[tuple[Any, Any]] = set()
    limit = len(graph)
    for edge in sorted(reticulate, key=lambda e: (depth[e[1]], depth[e[0]])):
        trial = horizontal | {edge}
        levels = dict(depth)
        for _ in range(limit + 1):
            changed = False
            for v in topo:
                for u in graph.predecessors(v):
                    need = levels[u] + (0 if (u, v) in trial else 1)
                    if levels[v] < need:
                        levels[v] = need
                        changed = True
            for p, h in trial:
                if levels[p] < levels[h]:
                    levels[p] = levels[h]
                    changed = True
            if not changed:
                break
            if max(levels.values()) > limit:
                break
        else:
            continue
        if changed:  # bound exceeded: positive cycle
            continue
        horizontal = trial
        depth.update(levels)
    return horizontal


def _preorder(root: Any, kids: Children) -> list[Any]:
    """Nodes of the backbone tree in preorder, children left to right."""
    order = []
    stack = [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(reversed(kids[n]))
    return order


def _x_coordinates(root: Any, kids: Children, leaf_gap: float) -> dict[Any, float]:
    """Leaves left to right in backbone order; internal nodes above the mean of their children."""
    order = _preorder(root, kids)
    x: dict[Any, float] = {}
    next_x = 0.0
    for n in order:
        if not kids[n]:
            x[n] = next_x
            next_x += leaf_gap
    for n in reversed(order):
        if kids[n]:
            x[n] = sum(x[c] for c in kids[n]) / len(kids[n])
    return x


def _barycenter_pass(
    root: Any, kids: Children, reticulate: list[tuple[Any, Any]], leaf_gap: float
) -> None:
    """Sort siblings by the mean x of the reticulate partners of their subtrees."""
    x = _x_coordinates(root, kids, leaf_gap)
    pull_sum = {n: 0.0 for n in kids}
    pull_count = {n: 0 for n in kids}
    for u, v in reticulate:
        pull_sum[u] += x[v]
        pull_count[u] += 1
        pull_sum[v] += x[u]
        pull_count[v] += 1
    for n in reversed(_preorder(root, kids)):
        for c in kids[n]:
            pull_sum[n] += pull_sum[c]
            pull_count[n] += pull_count[c]
    for cs in kids.values():
        if len(cs) > 1:
            cs.sort(key=lambda c: pull_sum[c] / pull_count[c] if pull_count[c] else x[c])


def _score(
    root: Any,
    kids: Children,
    edges: list[tuple[Any, Any]],
    reticulate: list[tuple[Any, Any]],
    depth: dict[Any, int],
    leaf_gap: float,
) -> tuple[int, float]:
    """(number of edge crossings, total horizontal span of reticulate edges); lower is better."""
    x = _x_coordinates(root, kids, leaf_gap)
    segments = np.array([(x[u], depth[u], x[v], depth[v]) for u, v in edges], dtype=float)
    span = sum(abs(x[u] - x[v]) for u, v in reticulate)
    return _count_crossings(segments), round(span, 9)


def _local_search(
    root: Any,
    kids: Children,
    edges: list[tuple[Any, Any]],
    reticulate: list[tuple[Any, Any]],
    depth: dict[Any, int],
    leaf_gap: float,
    score: Callable[[Children], tuple[int, float]] | None = None,
    nodes: list[Any] | None = None,
) -> tuple[int, float]:
    """Swap adjacent siblings (of ``nodes``, default all) while that improves the score; returns the final score."""
    if score is None:

        def score(k: Children) -> tuple[int, float]:
            return _score(root, k, edges, reticulate, depth, leaf_gap)

    best = score(kids)
    improved = True
    candidates = list(kids) if nodes is None else nodes
    while improved:
        improved = False
        for n in candidates:
            cs = kids[n]
            for i in range(len(cs) - 1):
                cs[i], cs[i + 1] = cs[i + 1], cs[i]
                candidate = score(kids)
                if candidate < best:
                    best = candidate
                    improved = True
                else:
                    cs[i], cs[i + 1] = cs[i + 1], cs[i]
    return best
