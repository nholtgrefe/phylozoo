"""
Shared layout utilities for PhyloZoo visualization.

This module provides common layout computation (NetworkX/Graphviz dispatch)
and position normalization used by all layout modules. Layouts determine
only node positions (placement); node sizes come from the style.
"""

from __future__ import annotations

from typing import Any, TypeVar

import networkx as nx
import numpy as np

from phylozoo.utils.exceptions import (
    PhyloZooImportError,
    PhyloZooLayoutError,
)

T = TypeVar("T")

GRAPHVIZ_LAYOUTS = ("dot", "neato", "fdp", "sfdp", "twopi", "circo")

NETWORKX_LAYOUTS = {
    "spring": nx.spring_layout,
    "circular": nx.circular_layout,
    "kamada_kawai": nx.kamada_kawai_layout,
    "planar": nx.planar_layout,
    "random": nx.random_layout,
    "shell": nx.shell_layout,
    "spectral": nx.spectral_layout,
    "spiral": nx.spiral_layout,
    "bipartite": nx.bipartite_layout,
}


def compute_nx_positions(
    G: nx.Graph,
    layout: str = "spring",
    **kwargs: Any,
) -> dict[Any, tuple[float, float]]:
    """
    Compute node positions using NetworkX or Graphviz.

    Parameters
    ----------
    G : nx.Graph
        NetworkX graph (DiGraph, Graph, MultiGraph, etc.).
    layout : str, optional
        Layout algorithm name. NetworkX: 'spring', 'circular', 'kamada_kawai',
        'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite'.
        Graphviz: 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'spring'.
    **kwargs
        Additional parameters passed to the layout algorithm.

    Returns
    -------
    dict
        Node ID -> (x, y) position mapping.

    Raises
    ------
    PhyloZooLayoutError
        If layout is not supported or computation fails.
    PhyloZooImportError
        If Graphviz layout requested but pygraphviz not installed.
    """
    pos: dict[Any, tuple[float, float]]

    if layout in GRAPHVIZ_LAYOUTS:
        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog=layout, **kwargs)
        except ImportError:
            raise PhyloZooImportError(
                f"Graphviz layout '{layout}' requires pygraphviz. "
                "Install with: pip install pygraphviz"
            )
        except Exception as e:
            raise PhyloZooLayoutError(f"Graphviz layout '{layout}' failed: {e}") from e
    elif layout in NETWORKX_LAYOUTS:
        pos = NETWORKX_LAYOUTS[layout](G, **kwargs)
    else:
        supported = ", ".join(sorted(NETWORKX_LAYOUTS.keys()) + list(GRAPHVIZ_LAYOUTS))
        raise PhyloZooLayoutError(
            f"Unsupported layout algorithm: '{layout}'. " f"Supported: {supported}"
        )

    return pos


def normalize_positions(
    pos: dict[T, tuple[float, float]],
) -> dict[T, tuple[float, float]]:
    """
    Center positions at origin and scale to unit range.

    Parameters
    ----------
    pos : dict
        Node ID -> (x, y) position mapping.

    Returns
    -------
    dict
        Normalized positions (centered, max dimension = 1).
    """
    if not pos:
        return pos

    xs = [x for x, _ in pos.values()]
    ys = [y for _, y in pos.values()]
    if not xs or not ys:
        return pos

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    width = max_x - min_x if max_x != min_x else 1.0
    height = max_y - min_y if max_y != min_y else 1.0
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    scale = 1.0 / max(width, height) if max(width, height) > 0 else 1.0

    return {node: ((x - center_x) * scale, (y - center_y) * scale) for node, (x, y) in pos.items()}


def compute_layout_center(
    pos: dict[T, tuple[float, float]],
) -> tuple[float, float]:
    """
    Compute the center of a layout from node positions.

    Parameters
    ----------
    pos : dict
        Node ID -> (x, y) position mapping.

    Returns
    -------
    tuple[float, float]
        (center_x, center_y) of the bounding box.
    """
    if not pos:
        return (0.0, 0.0)
    xs = [x for x, _ in pos.values()]
    ys = [y for _, y in pos.values()]
    return ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)


def stress_majorization(
    graph: nx.Graph,
    pos: dict[Any, tuple[float, float]],
    iterations: int,
    lengths: dict[tuple[Any, Any], float] | float = 1.0,
) -> dict[Any, tuple[float, float]]:
    """
    Relax a drawing by stress majorization (SMACOF), the algorithm behind ``neato``.

    The target distance between two nodes is their shortest-path distance in
    the graph, with edges weighted by ``lengths``; pairs are weighted by the
    inverse square of that distance, so local structure dominates. Starting
    from a good drawing this converges to a nearby local optimum, which is why
    the layouts built on it are deterministic.

    Parameters
    ----------
    graph : nx.Graph
        Simple undirected graph of the drawing (must be connected).
    pos : dict
        Starting positions for every node of ``graph``.
    iterations : int
        Number of majorization steps.
    lengths : dict[tuple, float] | float, optional
        Target length per edge ``(u, v)`` (as iterated by ``graph.edges``), or
        one length for all edges. By default 1.0.

    Returns
    -------
    dict
        Node ID -> (x, y) position mapping.

    Examples
    --------
    >>> import math
    >>> import networkx as nx
    >>> g = nx.path_graph(3)
    >>> pos = stress_majorization(g, {0: (0.0, 0.0), 1: (1.0, 0.3), 2: (2.0, 0.0)}, 20)
    >>> round(math.dist(pos[0], pos[2]), 1)
    2.0
    """
    if len(pos) < 2 or iterations <= 0:
        return dict(pos)
    nodes = list(graph)
    index = {v: i for i, v in enumerate(nodes)}
    for u, v in graph.edges:
        graph.edges[u, v]["len"] = lengths if isinstance(lengths, float) else lengths[(u, v)]
    n = len(nodes)
    target = np.zeros((n, n))
    for u, dists in nx.all_pairs_dijkstra_path_length(graph, weight="len"):
        for v, d in dists.items():
            target[index[u], index[v]] = d
    weights = np.zeros_like(target)
    off = ~np.eye(n, dtype=bool)
    weights[off] = target[off] ** -2
    laplacian = np.diag(weights.sum(axis=1)) - weights
    solve = np.linalg.pinv(laplacian)

    coords = np.array([pos[v] for v in nodes], dtype=float)
    for _ in range(iterations):
        dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
        dist[~off] = 1.0
        dist[dist < 1e-9] = 1e-9
        b_matrix = -weights * target / dist
        b_matrix[~off] = 0.0
        np.fill_diagonal(b_matrix, -b_matrix.sum(axis=1))
        coords = solve @ (b_matrix @ coords)
    return {v: (float(coords[i, 0]), float(coords[i, 1])) for v, i in index.items()}


def count_crossings(segments: np.ndarray) -> int:
    """
    Number of properly crossing pairs among line segments.

    Parameters
    ----------
    segments : np.ndarray
        Array of shape (m, 4) with rows ``(x1, y1, x2, y2)``.

    Returns
    -------
    int
        Number of pairs of segments that cross at interior points (shared
        endpoints and collinear overlaps do not count).

    Examples
    --------
    >>> import numpy as np
    >>> count_crossings(np.array([[0, 0, 1, 1], [0, 1, 1, 0]]))
    1
    >>> count_crossings(np.array([[0, 0, 1, 0], [0, 1, 1, 1]]))
    0
    """
    m = len(segments)
    if m < 2:
        return 0
    i, j = np.triu_indices(m, 1)
    a, b = segments[i, :2], segments[i, 2:]
    c, d = segments[j, :2], segments[j, 2:]

    def orient(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> np.ndarray:
        cross = (q[:, 0] - p[:, 0]) * (r[:, 1] - p[:, 1]) - (q[:, 1] - p[:, 1]) * (
            r[:, 0] - p[:, 0]
        )
        return np.asarray(np.sign(cross))

    crossing = (orient(a, b, c) * orient(a, b, d) < 0) & (orient(c, d, a) * orient(c, d, b) < 0)
    return int(crossing.sum())
