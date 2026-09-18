"""
Tree-of-blobs placement shared by the ``pz-unrooted`` layouts of both network classes.

The network is first decomposed into its blobs (maximal cut-edge-free
subgraphs). Contracting every blob gives a tree, which is drawn with the
equal-angle algorithm for unrooted trees; each non-trivial blob is then drawn
separately as a regular polygon inside the space reserved for it. The order of
the pendant subtrees around a blob follows the cyclic order of their
attachment nodes on the blob's outer face, so cut edges never cross blob
edges. Two optional finishing steps spread the drawing: equal-daylight rounds
(rigid rotations, exact edge lengths) and stress majorization.
"""

from __future__ import annotations

import math
from typing import Any, Callable

import networkx as nx
import numpy as np

from phylozoo.utils.exceptions import PhyloZooLayoutError
from phylozoo.viz._layout_utils import count_crossings, normalize_positions, stress_majorization

Point = tuple[float, float]

_TWO_PI = 2.0 * math.pi
# Angular space claimed by a blob node on its circle, in "leaf" units.
_NODE_WEIGHT = 1.0


def simple_graph(network: Any) -> nx.Graph:
    """
    Simple undirected graph of a network: directions dropped, parallel edges counted in ``mult``.

    Parameters
    ----------
    network : DirectedPhyNetwork | SemiDirectedPhyNetwork
        Any network whose ``_graph`` exposes ``nodes`` and ``edges``.

    Returns
    -------
    nx.Graph
        The underlying simple graph.
    """
    graph: nx.Graph = nx.Graph()
    graph.add_nodes_from(network._graph.nodes)
    for u, v in network._graph.edges:
        if u == v:
            continue
        if graph.has_edge(u, v):
            graph.edges[u, v]["mult"] += 1
        else:
            graph.add_edge(u, v, mult=1)
    return graph


def tree_of_blobs_positions(
    graph: nx.Graph,
    leaves: set[Any],
    edge_length: float = 1.0,
    node_spacing: float = 0.7,
    inner_scale: float = 0.6,
    start_angle: float = 0.0,
    refine: int = 50,
    daylight: int = 0,
) -> dict[Any, Point]:
    """
    Place the nodes of a connected graph with the tree-of-blobs method.

    See :func:`phylozoo.viz.sdnetwork.layout.unrooted.compute_pz_unrooted_layout`
    for the algorithm and the meaning of the parameters.

    Parameters
    ----------
    graph : nx.Graph
        Simple undirected graph (as from :func:`simple_graph`).
    leaves : set
        Leaf nodes (they weight the angular wedges).
    edge_length, node_spacing, inner_scale, start_angle, refine, daylight
        Layout parameters.

    Returns
    -------
    dict
        Node ID -> (x, y), centred and scaled to unit range.

    Raises
    ------
    PhyloZooLayoutError
        If the graph is empty or not connected.
    """
    if len(graph) == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")
    if not nx.is_connected(graph):
        raise PhyloZooLayoutError("pz-unrooted layout requires a connected network")
    # Cut edges (an edge with parallel copies is never a cut edge) and blobs.
    bridges = {(u, v) for u, v in nx.bridges(graph) if graph.edges[u, v].get("mult", 1) == 1}
    without = graph.copy()
    without.remove_edges_from(bridges)
    blob_sets: list[set[Any]] = [set(c) for c in nx.connected_components(without)]
    blob_of = {v: i for i, b in enumerate(blob_sets) for v in b}
    tree: nx.Graph = nx.Graph()
    tree.add_nodes_from(range(len(blob_sets)))
    for u, v in bridges:
        bu, bv = blob_of[u], blob_of[v]
        tree.add_edge(bu, bv, anchors={bu: u, bv: v})

    # Space demanded by each blob: its leaves, plus one leaf-equivalent per node
    # of a non-trivial blob (those nodes need room on the blob's circle).
    demand = [len(b & leaves) + (_NODE_WEIGHT * len(b) if len(b) > 1 else 0.0) for b in blob_sets]
    root, parent, children, weight = _root_at_centroid(tree, demand, blob_sets)

    def anchor(b: int, c: int) -> Any:
        """Node of blob ``b`` incident to the cut edge towards blob ``c``."""
        return tree.edges[b, c]["anchors"][b]

    positions: dict[Any, Point] = {}

    def fan_out(kids: list[int], origin: Point, lo: float, hi: float, scale: float) -> None:
        """Give each child a sub-wedge of [lo, hi] proportional to its weight."""
        total = sum(weight[c] for c in kids)
        a = lo
        for c in kids:
            width = (hi - lo) * weight[c] / total
            place(c, a, a + width, origin, scale)
            a += width

    def place(b: int, lo: float, hi: float, origin: Point | None, scale: float) -> None:
        """Place blob ``b`` and, recursively, its subtree inside the wedge [lo, hi]."""
        nodes = blob_sets[b]
        theta = (lo + hi) / 2
        length = edge_length * scale
        if len(nodes) == 1:
            center = _polar(origin, length, theta)
            positions[next(iter(nodes))] = center
            if children[b]:
                fan_out(children[b], center, lo, hi, scale)
            return

        blob_graph = nx.Graph(graph.subgraph(nodes))
        parent_blob = parent[b]
        parent_anchor = anchor(b, parent_blob) if parent_blob is not None else None
        kids_at: dict[Any, list[int]] = {}
        for c in children[b]:
            kids_at.setdefault(anchor(b, c), []).append(c)
        ring = _outer_cycle(
            blob_graph, {a: sum(weight[c] for c in cs) for a, cs in kids_at.items()}, parent_anchor
        )
        if parent_anchor is not None:
            i = ring.index(parent_anchor)
            ring = ring[i + 1 :] + ring[:i]

        # Angular slots on the circle: the parent anchor faces the parent, the
        # other outer nodes are evenly spaced (an evenly sided polygon); subtrees
        # on the parent anchor itself get a share at the edge of the span.
        anchor_kids = kids_at.get(parent_anchor, []) if parent_anchor is not None else []
        # The ring is a regular polygon (one slot per node, plus one for the parent
        # anchor); the wedge only orients it, unless the wedge is wider still.
        span = max(
            hi - lo, _TWO_PI * len(ring) / (len(ring) + (1 if parent_anchor is not None else 0))
        )
        ring_lo, ring_hi = theta - span / 2, theta + span / 2
        share = (ring_hi - ring_lo) * sum(weight[c] for c in anchor_kids)
        share /= sum(weight[c] for c in anchor_kids) + len(ring)
        angles: dict[Any, float] = {}
        slots: dict[Any, tuple[float, float]] = {}
        a = ring_lo + share
        for v in ring:
            width = (ring_hi - ring_lo - share) / len(ring)
            slots[v] = (a, a + width)
            angles[v] = a + width / 2
            a += width
        if parent_anchor is not None:
            angles[parent_anchor] = theta + math.pi

        radius = _blob_radius(list(angles.values()), len(nodes) - len(angles), node_spacing)
        center = _polar(origin, length + radius, theta)
        for v, ang in angles.items():
            positions[v] = _polar(center, radius, ang)
        interior = [v for v in nodes if v not in angles]
        if interior:
            circle = {v: positions[v] for v in angles}
            positions.update(_tutte_positions(blob_graph, circle, interior, node_spacing))

        for v in ring:
            if v in kids_at:
                fan_out(kids_at[v], positions[v], *slots[v], scale)
        if anchor_kids:
            # Fan them out just beside the anchor's blob edge to ring[0], between
            # that edge and the edge to the parent: outside the convex polygon,
            # and as close to their reserved share of the wedge as possible.
            px, py = positions[parent_anchor]
            beta = math.atan2(positions[ring[0]][1] - py, positions[ring[0]][0] - px)
            gap = (beta - (theta + math.pi)) % _TWO_PI
            margin = min(0.3, gap / 4)
            width = min(share, gap - 2 * margin)
            fan_out(anchor_kids, (px, py), beta - margin - width, beta - margin, scale)
        for v, kids in kids_at.items():
            if v in angles:
                continue
            # Interior anchor: point the subtree into the largest gap between the
            # node's incident edges (the widest face around it).
            px, py = positions[v]
            incident = sorted(
                math.atan2(positions[u][1] - py, positions[u][0] - px)
                for u in blob_graph.neighbors(v)
            )
            gaps = [
                (incident[(i + 1) % len(incident)] - incident[i]) % _TWO_PI or _TWO_PI
                for i in range(len(incident))
            ]
            i = max(range(len(gaps)), key=gaps.__getitem__)
            lo_gap, width = incident[i], gaps[i]
            fan_out(kids, (px, py), lo_gap + 0.2 * width, lo_gap + 0.8 * width, scale * inner_scale)

    place(root, start_angle, start_angle + _TWO_PI, None, 1.0)
    if daylight > 0:
        movable = {frozenset(e) for e in bridges}
        rigid = [b for b in blob_sets if len(b) > 1]
        positions = _guarded(
            graph,
            positions,
            daylight,
            lambda pos, n: _equal_daylight(graph, pos, movable, rigid, n),
            1,
        )
    if refine > 0 and len(positions) > 2:
        # Blob edges target node_spacing (evenly sided polygons), cut edges their drawn length.
        lengths = {
            (u, v): (
                node_spacing if blob_of[u] == blob_of[v] else math.dist(positions[u], positions[v])
            )
            for u, v in graph.edges
        }
        positions = _guarded(
            graph,
            positions,
            refine,
            lambda pos, n: stress_majorization(graph, pos, n, lengths),
            _STRESS_CHUNK,
        )
    return normalize_positions(positions)


def _polar(origin: Point | None, distance: float, angle: float) -> Point:
    """Point at ``distance`` from ``origin`` in direction ``angle`` (origin None = (0, 0))."""
    if origin is None:
        return (0.0, 0.0)
    return (origin[0] + distance * math.cos(angle), origin[1] + distance * math.sin(angle))


def _root_at_centroid(
    tree: nx.Graph, demand: list[float], blob_sets: list[set[Any]]
) -> tuple[int, dict[int, int | None], dict[int, list[int]], list[float]]:
    """
    Root the tree of blobs at its centroid.

    Returns the root, the parent map, the ordered children lists and the total
    demand of every subtree (at least 1, so every subtree gets some angle).
    """

    def rooted(root: int) -> tuple[dict[int, int | None], list[int], list[float]]:
        parent: dict[int, int | None] = {root: None}
        order = [root]
        for u, v in nx.bfs_edges(tree, root):
            parent[v] = u
            order.append(v)
        sub = list(demand)
        for v in reversed(order):
            p = parent[v]
            if p is not None:
                sub[p] += sub[v]
        return parent, order, sub

    parent, order, sub = rooted(next(iter(tree.nodes)))
    total = sum(demand)

    def heaviest_side(b: int) -> float:
        parts = [sub[c] for c in tree.neighbors(b) if parent[c] == b]
        if parent[b] is not None:
            parts.append(total - sub[b])
        return max(parts, default=0.0)

    root = min(tree.nodes, key=lambda b: (heaviest_side(b), -len(blob_sets[b]), b))
    parent, order, sub = rooted(root)
    children: dict[int, list[int]] = {b: [] for b in tree.nodes}
    for v in order[1:]:
        children[parent[v]].append(v)  # type: ignore[index]
    return root, parent, children, [max(w, 1.0) for w in sub]


def _outer_cycle(graph: nx.Graph, anchor_weight: dict[Any, float], parent_anchor: Any) -> list[Any]:
    """
    Cyclic order of the nodes drawn on the blob's circle.

    For a planar blob this is the face of a planar embedding that contains the
    parent anchor and the heaviest pendant subtrees. Otherwise all anchors are
    used, ordered by their angle in a force-directed drawing of the blob.
    """
    is_planar, embedding = nx.check_planarity(graph)
    if is_planar:
        seen: set[tuple[Any, Any]] = set()
        faces: list[list[Any]] = []
        for u, v in embedding.edges():
            if (u, v) not in seen:
                walk = embedding.traverse_face(u, v, mark_half_edges=seen)
                faces.append(list(dict.fromkeys(walk)))

        def score(face: list[Any]) -> tuple[bool, float, int]:
            return (
                parent_anchor in face,
                sum(anchor_weight.get(v, 0.0) for v in face),
                len(face),
            )

        return max(faces, key=score)

    pos = nx.spring_layout(graph, seed=0)
    cx, cy = np.mean(list(pos.values()), axis=0)
    outer = [v for v in graph if v in anchor_weight or v == parent_anchor]
    return sorted(outer, key=lambda v: math.atan2(pos[v][1] - cy, pos[v][0] - cx))


def _blob_radius(angles: list[float], n_interior: int, spacing: float) -> float:
    """Circle radius that keeps consecutive outer nodes about ``spacing`` apart."""
    n = len(angles)
    regular = spacing / (2 * math.sin(math.pi / n)) if n > 1 else 0.0
    ordered = sorted(a % _TWO_PI for a in angles)
    gaps = [(ordered[(i + 1) % n] - ordered[i]) % _TWO_PI for i in range(n)]
    min_gap = min(g for g in gaps if g > 1e-9) if n > 1 else math.pi
    crowded = spacing / (2 * math.sin(min(min_gap, math.pi) / 2))
    radius = min(max(regular, crowded), 3 * regular)
    if n_interior:
        radius = max(radius, spacing * math.sqrt(n_interior + 1))
    return radius


def _tutte_positions(
    graph: nx.Graph, fixed: dict[Any, Point], interior: list[Any], spacing: float
) -> dict[Any, Point]:
    """
    Place interior nodes at the barycentre of their neighbours (Tutte embedding).

    If that collapses nodes onto each other (the blob is not 3-connected), a
    short force-directed refinement with the circle nodes fixed is applied.
    """
    index = {v: i for i, v in enumerate(interior)}
    laplacian = np.zeros((len(interior), len(interior)))
    rhs = np.zeros((len(interior), 2))
    for v, i in index.items():
        for u in graph.neighbors(v):
            laplacian[i, i] += 1
            if u in index:
                laplacian[i, index[u]] -= 1
            else:
                rhs[i] += fixed[u]
    solution = np.linalg.solve(laplacian, rhs)
    pos = {v: (float(solution[i, 0]), float(solution[i, 1])) for v, i in index.items()}

    everything = dict(fixed)
    everything.update(pos)
    points = np.array([everything[v] for v in graph])
    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    np.fill_diagonal(distances, np.inf)
    if distances.min() < 0.3 * spacing:
        refined = nx.spring_layout(
            graph, pos=everything, fixed=list(fixed), k=spacing, iterations=50, seed=0
        )
        pos = {v: (float(refined[v][0]), float(refined[v][1])) for v in interior}
    return pos


def _equal_daylight(
    graph: nx.Graph,
    pos: dict[Any, Point],
    movable: set[frozenset[Any]],
    blobs_: list[set[Any]],
    rounds: int,
) -> dict[Any, Point]:
    """
    Equal-daylight adjustment: rotate pendant subtrees so the free angles are equal.

    Every pivot is either a tree node (all incident edges are cut edges) or a
    whole non-trivial blob. The components of the graph without the pivot are
    rigid subtrees (blobs included) hanging off it by one cut edge each. Their
    angular extents are measured from the pivot's centre and the subtrees are
    rotated about their attachment node so that the empty angle between
    consecutive subtrees is the same all around. For a blob, a subtree is never
    rotated more than 80 degrees away from the outward normal at its
    attachment node, so its cut edge stays outside the polygon.
    """
    pos = dict(pos)
    pivots: list[set[Any]] = [set(b) for b in blobs_]
    pivots += [
        {v}
        for v in graph
        if graph.degree(v) >= 2 and all(frozenset((u, v)) in movable for u in graph.neighbors(v))
    ]
    for _ in range(rounds):
        for pivot in pivots:
            cx = sum(pos[v][0] for v in pivot) / len(pivot)
            cy = sum(pos[v][1] for v in pivot) / len(pivot)
            rest = graph.subgraph(n for n in graph if n not in pivot)
            parts = []  # (anchor angle, extent lo, extent hi, nodes, anchor)
            for a in pivot:
                for u in graph.neighbors(a):
                    if u in pivot:
                        continue
                    nodes = nx.node_connected_component(rest, u)
                    direction = math.atan2(pos[u][1] - cy, pos[u][0] - cx)
                    rel = [
                        (
                            (math.atan2(pos[n][1] - cy, pos[n][0] - cx) - direction + math.pi)
                            % _TWO_PI
                        )
                        - math.pi
                        for n in nodes
                    ]
                    parts.append((direction, min(rel), max(rel), nodes, a))
            parts.sort(key=lambda p: p[0])
            if len(parts) < 2:
                continue
            gap = (_TWO_PI - sum(hi - lo for _, lo, hi, _, _ in parts)) / len(parts)
            if gap <= 0:
                continue
            cursor = parts[0][0] + parts[0][1]
            for direction, lo, hi, nodes, a in parts:
                shift = cursor - (direction + lo)
                ax_, ay_ = pos[a]
                if len(pivot) > 1:
                    # Keep the cut edge within 80 degrees of the outward normal at the anchor.
                    normal = math.atan2(ay_ - cy, ax_ - cx)
                    first = next(n for n in nodes if graph.has_edge(a, n))
                    edge = math.atan2(pos[first][1] - ay_, pos[first][0] - ax_)
                    off = ((edge + shift - normal + math.pi) % _TWO_PI) - math.pi
                    limit = math.radians(80)
                    shift += max(-limit, min(limit, off)) - off
                c, s = math.cos(shift), math.sin(shift)
                for n in nodes:
                    dx, dy = pos[n][0] - ax_, pos[n][1] - ay_
                    pos[n] = (ax_ + c * dx - s * dy, ay_ + s * dx + c * dy)
                cursor += (hi - lo) + gap
    return pos


_STRESS_CHUNK = 10


def _guarded(
    graph: nx.Graph,
    pos: dict[Any, Point],
    total: int,
    step: Callable[[dict[Any, Point], int], dict[Any, Point]],
    chunk: int,
) -> dict[Any, Point]:
    """
    Apply an iterative relaxation without ever increasing the number of edge crossings.

    ``step(pos, n)`` performs ``n`` iterations from ``pos``. The iterations run
    in chunks; a drawing is kept only if it has at most as many crossings as
    the starting drawing, otherwise the iteration continues from it but the
    last accepted drawing is what is returned. Neither stress majorization nor
    the daylight rotations know about planarity, so without this guard they
    can fold a crossing-free tree-of-blobs drawing.
    """
    edges = list(graph.edges)

    def crossings(p: dict[Any, Point]) -> int:
        return count_crossings(np.array([(*p[u], *p[v]) for u, v in edges], dtype=float))

    limit = crossings(pos)
    best, current = pos, pos
    for start in range(0, total, chunk):
        current = step(current, min(chunk, total - start))
        if crossings(current) <= limit:
            best = current
    return best
