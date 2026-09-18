"""
Pretty-printing: text drawings of phylogenetic networks.

A network is drawn root-left / leaves-right on a character grid, using the
rectangular ``pz-cladogram`` placement: tree edges as box-drawing lines, hybrid
edges as dotted lines with an arrowhead at the hybrid node. Semi-directed
networks are rooted first. Needs no matplotlib.

Legend::

    ○  root          ●  tree node / leaf      ◆  hybrid node
    ─ │ ┌ └ ├ …  tree edges     ┄ ┊  hybrid edges, ending in > < v ^
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

from phylozoo.utils.exceptions import (
    PhyloZooEmptyNetworkWarning,
    PhyloZooSingleNodeNetworkWarning,
    PhyloZooTypeError,
    PhyloZooValueError,
)

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

_BOX = {
    frozenset("EW"): "─",
    frozenset("NS"): "│",
    frozenset("SE"): "┌",
    frozenset("SW"): "┐",
    frozenset("NE"): "└",
    frozenset("NW"): "┘",
    frozenset("NSE"): "├",
    frozenset("NSW"): "┤",
    frozenset("SEW"): "┬",
    frozenset("NEW"): "┴",
    frozenset("NSEW"): "┼",
    frozenset("E"): "─",
    frozenset("W"): "─",
    frozenset("N"): "│",
    frozenset("S"): "│",
}
_DOTTED_H, _DOTTED_V, _CROSS = "┄", "┊", "┼"
_ROOT, _NODE, _HYBRID, _LEAF = "○", "●", "◆", "●"
_HEADS = {"E": ">", "W": "<", "S": "v", "N": "^"}


def to_pretty_print(
    network: "DirectedPhyNetwork | SemiDirectedPhyNetwork",
    col_width: int | None = None,
    rows_per_leaf: int | None = None,
    max_width: int = 120,
    max_leaves: int = 80,
    root_location: Any = None,
    **layout_kwargs: Any,
) -> str:
    """
    Draw a network as text, root on the left and leaves on the right.

    Tree edges are box-drawing lines (``─ │ ┌ └ ├`` ...), hybrid edges dotted
    lines (``┄ ┊``) ending in an arrowhead at the hybrid node (``◆``); the
    root is ``○``, other nodes ``●``, and leaf labels follow their leaf. The
    placement is that of the rectangular ``pz-cladogram`` layout, so hybrid edges
    are horizontal or vertical where the network allows it and sibling
    subtrees stay together. A semi-directed network is rooted first with
    :func:`~phylozoo.core.network.sdnetwork.derivations.to_d_network`.

    Parallel edges are drawn as a loop: the second copy detours one row above
    (or below) the first, so both stay visible.

    Large networks are compacted automatically: the column width shrinks
    until the drawing fits ``max_width`` characters (down to 3) and the blank
    row between leaves is dropped above 25 leaves. Networks with more than
    ``max_leaves`` leaves are refused, because one text row per leaf stops
    being readable; raise ``max_leaves`` or use :func:`~phylozoo.viz.plot`.

    Parameters
    ----------
    network : DirectedPhyNetwork | SemiDirectedPhyNetwork
        The network to draw.
    col_width : int | None, optional
        Characters per layer. None picks the largest of 8, 6, 4, 3 that fits
        ``max_width``. By default None.
    rows_per_leaf : int | None, optional
        Rows per leaf: 2 leaves a blank row between leaves, 1 is compact.
        None chooses 2 up to 25 leaves and 1 above. By default None.
    max_width : int, optional
        Target maximum line length used to pick ``col_width``. By default 120.
    max_leaves : int, optional
        Refuse networks with more leaves than this. By default 80.
    root_location : node or edge, optional
        Semi-directed networks only: where to root (see ``to_d_network``).
        None lets ``to_d_network`` choose. By default None.
    **layout_kwargs
        Options of :func:`~phylozoo.viz.dnetwork.layout.cladogram.compute_pz_cladogram_layout`
        (e.g. ``trials``, ``seed``, ``horizontal_reticulations``).

    Returns
    -------
    str
        The drawing, lines separated by newlines, without trailing spaces.

    Raises
    ------
    PhyloZooTypeError
        If ``network`` is not a directed or semi-directed network.
    PhyloZooValueError
        If the network has more than ``max_leaves`` leaves.

    Examples
    --------
    >>> from phylozoo import DirectedPhyNetwork
    >>> from phylozoo.viz import to_pretty_print
    >>> net = DirectedPhyNetwork(
    ...     edges=[(5, 3), (5, 4), (3, 2), (4, 2), (2, 1), (3, 6), (4, 7)],
    ...     nodes=[(1, {'label': 'A'}), (6, {'label': 'B'}), (7, {'label': 'C'})]
    ... )
    >>> print(to_pretty_print(net))
            ┌───────────────● B
    ┌───────●
    ○       └┄┄┄┄┄┄>◆───────● A
    │               ^
    └───────────────●───────● C
    """
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

    if isinstance(network, DirectedPhyNetwork):
        rooted = network
    elif isinstance(network, SemiDirectedPhyNetwork):
        from phylozoo.core.network.sdnetwork.derivations import to_d_network

        if network.number_of_nodes() <= 1:
            return _tiny(network)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=PhyloZooEmptyNetworkWarning)
            warnings.filterwarnings("ignore", category=PhyloZooSingleNodeNetworkWarning)
            rooted = to_d_network(network, root_location=root_location)
    else:
        raise PhyloZooTypeError(
            f"to_pretty_print expects a DirectedPhyNetwork or SemiDirectedPhyNetwork, got {type(network).__name__}"
        )
    if rooted.number_of_nodes() <= 1:
        return _tiny(rooted)
    n_leaves = len(rooted.leaves)
    if n_leaves > max_leaves:
        raise PhyloZooValueError(
            f"Network has {n_leaves} leaves, more than max_leaves={max_leaves}; a text drawing "
            "needs one row per leaf. Raise max_leaves or use phylozoo.viz.plot()."
        )
    return _draw(rooted, col_width, rows_per_leaf, max_width, layout_kwargs)


def _tiny(network: Any) -> str:
    """Drawing of an empty or single-node network."""
    nodes = list(network._graph.nodes)
    if not nodes:
        return ""
    return f"{_LEAF} {network.get_label(nodes[0]) or nodes[0]}"


def _draw(
    net: "DirectedPhyNetwork",
    col_width: int | None,
    rows_per_leaf: int | None,
    max_width: int,
    layout_kwargs: dict[str, Any],
) -> str:
    from phylozoo.viz.dnetwork.layout.cladogram import compute_pz_cladogram_layout

    layout_kwargs.setdefault("rectangular", True)
    layout = compute_pz_cladogram_layout(net, direction="LR", **layout_kwargs)
    pos = layout.positions
    layers = sorted({round(x, 9) for x, _ in pos.values()})
    leaves = sorted((n for n in pos if n in net.leaves), key=lambda n: -pos[n][1])
    labels = [str(net.get_label(n) or n) for n in leaves]
    if rows_per_leaf is None:
        rows_per_leaf = 2 if len(leaves) <= 25 else 1
    if col_width is None:
        longest = max(len(s) for s in labels) + 2
        col_width = next(
            (w for w in (8, 6, 4, 3) if (len(layers) - 1) * w + 1 + longest <= max_width), 3
        )

    top = pos[leaves[0]][1]
    step = (top - pos[leaves[-1]][1]) / max(len(leaves) - 1, 1) or 1.0
    col_of = {x: i * col_width for i, x in enumerate(layers)}

    def cell(point: tuple[float, float]) -> tuple[int, int]:
        x, y = point
        nearest = min(layers, key=lambda lx: abs(lx - x))
        return col_of[nearest], round((top - y) / step * rows_per_leaf)

    n_rows = max(cell(p)[1] for p in pos.values()) + 1
    n_cols = max(cell(p)[0] for p in pos.values()) + 1
    tree: list[list[set[str]]] = [[set() for _ in range(n_cols)] for _ in range(n_rows)]
    dotted: list[list[set[str]]] = [[set() for _ in range(n_cols)] for _ in range(n_rows)]

    def segment(grid: list[list[set[str]]], a: tuple[int, int], b: tuple[int, int]) -> None:
        (c0, r0), (c1, r1) = a, b
        if r0 == r1:
            for c in range(min(c0, c1), max(c0, c1) + 1):
                if c > min(c0, c1):
                    grid[r0][c].add("W")
                if c < max(c0, c1):
                    grid[r0][c].add("E")
        else:
            for r in range(min(r0, r1), max(r0, r1) + 1):
                if r > min(r0, r1):
                    grid[r][c0].add("N")
                if r < max(r0, r1):
                    grid[r][c0].add("S")

    heads: list[tuple[tuple[int, int], tuple[int, int]]] = []
    copies: dict[tuple[Any, Any], list[int]] = {}
    for u, v, key in layout.edge_routes:
        copies.setdefault((u, v), []).append(key)
    for (u, v, key), route in layout.edge_routes.items():
        cells = [cell(p) for p in route.points]
        cells = [c for i, c in enumerate(cells) if i == 0 or c != cells[i - 1]]
        nth = sorted(copies[(u, v)]).index(key)
        if nth:  # extra parallel copy: detour one row per copy so it stays visible
            shift = -nth if cells[0][1] - nth >= 0 else nth
            cells = (
                [cells[0], (cells[0][0], cells[0][1] + shift)]
                + [(c, r + shift) for c, r in cells[1:-1]]
                + [(cells[-1][0], cells[-1][1] + shift), cells[-1]]
            )
            cells = [c for i, c in enumerate(cells) if i == 0 or c != cells[i - 1]]
        grid = dotted if route.edge_type.is_hybrid else tree
        for a, b in zip(cells, cells[1:]):
            segment(grid, a, b)
        if route.edge_type.is_hybrid and len(cells) > 1:
            heads.append((cells[-1], cells[-2]))

    out = [[" "] * n_cols for _ in range(n_rows)]
    for r in range(n_rows):
        for c in range(n_cols):
            solid, dots = tree[r][c], dotted[r][c]
            if solid and dots:
                # Same orientation: the hybrid edge runs along the tree edge, show it dotted;
                # different orientation: a real crossing.
                parallel = bool(solid & {"N", "S"}) == bool(dots & {"N", "S"})
                out[r][c] = (_DOTTED_V if dots & {"N", "S"} else _DOTTED_H) if parallel else _CROSS
            elif solid:
                out[r][c] = _BOX[frozenset(solid)]
            elif dots & {"N", "S"} and dots & {"E", "W"}:
                out[r][c] = _BOX[frozenset(dots)]  # corner of a dotted detour
            elif dots:
                out[r][c] = _DOTTED_V if dots & {"N", "S"} else _DOTTED_H
    for (c, r), (pc, pr) in heads:  # arrowhead in the cell before the hybrid node
        if pr == r:
            direction = "E" if pc < c else "W"
            ac, ar = (c - 1, r) if pc < c else (c + 1, r)
        else:
            direction = "S" if pr < r else "N"
            ac, ar = (c, r - 1) if pr < r else (c, r + 1)
        if 0 <= ar < n_rows and 0 <= ac < n_cols and not tree[ar][ac]:
            out[ar][ac] = _HEADS[direction]  # never over a tree edge
    for node, point in pos.items():
        c, r = cell(point)
        if node in net.leaves:
            out[r][c] = _LEAF
        elif node in net.hybrid_nodes:
            out[r][c] = _HYBRID
        elif node == net.root_node:
            out[r][c] = _ROOT
        else:
            out[r][c] = _NODE
    for leaf, label in zip(leaves, labels):
        out[cell(pos[leaf])[1]] += list(" " + label)
    return "\n".join("".join(row).rstrip() for row in out)
